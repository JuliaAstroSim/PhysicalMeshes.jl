# Halo exchange for `DistributedMesh`, built on `ParallelOperations`.
#
# Each worker has a module-level `HaloStash.incoming` dict keyed by
# `(sender_rank::Int, field_name::Symbol)`. Halo exchange has three
# phases:
#
#   1. send     every process fires `ParallelOperations.sendto_async`
#               to deposit its boundary slice into each neighbour's
#               `HaloStash.incoming`. The future is collected so we
#               can later wait for it.
#   2. barrier  `ParallelOperations.gather` over (self ∪ neighbours).
#               Once it returns on every participant, every deposit
#               from phase 1 is guaranteed visible on the receiver.
#   3. apply    every process reads its own `HaloStash.incoming` and
#               broadcasts the slice into the local ghost zone via
#               `selectdim`.
#
# `sendto_async` ultimately runs `Core.eval(:(LHS = data))` on the
# remote worker via `@spawnat`, so the assignment happens in the
# *receiver* module and is process-global — exactly the property the
# old `TaskLocalStorage` version lacked.
#
# Direction convention
# --------------------
# `direction ∈ {±1, ±2, ±3}` encodes axis `abs(direction)` and side
# (`+` = max-end / high index, `-` = min-end / low index).
# `extract_halo(mesh, dir, field)` returns the slice at the side
# pointed to by `dir`; `apply_halo!(mesh, halo, dir, field)` writes
# the received slice into the *opposite* side of the ghost zone of
# `mesh`, since what is "min" on the neighbour is "max" from us.
# ---------------------------------------------------------------------------

using ParallelOperations: sendto_async, gather
using Distributed: Future

# Per-worker module-level stash for incoming halos.  Keyed by
# `(sender_rank::Int, field_name::Symbol)`, valued by the boundary
# slice as a regular `Array`.  Cleaned up after each apply to avoid
# unbounded growth across time steps.
module HaloStash
    const incoming = Dict{Tuple{Int,Symbol}, Any}()
end

# --- field-data adapter --------------------------------------------------
# Mirror the helper in `mesh/tools.jl` so this file is self-contained
# but stays consistent with the canonical extraction used elsewhere.
_field_data(f) = f
_field_data(f::ArrayScalarField) = f.data
_field_data(f::ArrayVectorField) = f.data

# --- direction -----------------------------------------------------------

"""
    _direction_to(rank, neighbor, size, dim) -> Int

Return the direction (axis + side) from `rank` to `neighbor` in a
partition of size `size` for a `dim`-dimensional grid. 1D and 3D use
the linear row layout; 2D uses the same `ceil(sqrt(size))` row /
`ceil(size/rows)` column layout as `partition_neighbors`.
"""
function _direction_to(rank::Int, neighbor::Int, size::Int, dim::Int)
    dim == 1 && return neighbor < rank ? -1 : 1
    dim == 3 && return neighbor < rank ? -1 : 1
    # 2D
    rows = Int(ceil(sqrt(size)))
    cols = Int(ceil(size / rows))
    rank_row, rank_col = divrem(rank, cols)
    neigh_row, neigh_col = divrem(neighbor, cols)
    if     neigh_row < rank_row; return -2   # up
    elseif neigh_row > rank_row; return  2   # down
    elseif neigh_col < rank_col; return -1   # left
    else;                        return  1   # right
    end
end

# --- extract / apply -----------------------------------------------------

"""
    _extract_halo(mesh, direction, field_name) -> AbstractArray

Return a copy of the boundary slice of `field_name` on the side of
`direction`. Uses `selectdim` to identify the slice but **materialises
to a fresh `Array`** (`collect`) before returning, so that the halo
is a plain owned array rather than a `SubArray` view. This matters
because the halo is later captured into a closure and shipped across
processes via `sendto_async` / `@spawnat`; `SubArray`s do not always
survive that round-trip cleanly (Julia 1.12 in particular can drop
the index information when serialising a `SubArray` embedded in an
`Expr`), which silently turns the halo into the wrong value on the
receiver. Owning the array up-front removes that ambiguity.
"""
function _extract_halo(mesh::DistributedMesh, direction::Int, field_name::Symbol)
    field = getproperty(mesh.local_mesh, field_name)
    data  = _field_data(field)
    NG    = mesh.local_mesh.config.NG
    ax    = abs(direction)
    n     = size(data, ax)
    slice = direction > 0 ? selectdim(data, ax, (n - NG + 1):n) :
                            selectdim(data, ax, 1:NG)
    return collect(slice)
end

"""
    _apply_halo!(mesh, halo, direction, field_name) -> mesh

Write the received `halo` into the ghost zone on the side pointed to
by `direction` (note: the *opposite* side of the neighbour's slice,
since min on the neighbour == max from us and vice versa).
"""
function _apply_halo!(mesh::DistributedMesh, halo::AbstractArray,
                      direction::Int, field_name::Symbol)
    field = getproperty(mesh.local_mesh, field_name)
    data  = _field_data(field)
    NG    = mesh.local_mesh.config.NG
    ax    = abs(direction)
    @assert size(halo, ax) == NG "halo thickness along axis $ax " *
                                "($(size(halo, ax))) ≠ NG=$NG"
    ghost = direction > 0 ? selectdim(data, ax, (size(data, ax) - NG + 1):size(data, ax)) :
                            selectdim(data, ax, 1:NG)
    ghost .= halo
    return mesh
end

# --- public API ----------------------------------------------------------

"""
$(TYPEDSIGNATURES)
Exchange boundary data for `field_name` (default `:rho`) with every
configured neighbour of `mesh`. No-op when the mesh has no neighbours
(serial run with a single rank).

The keyword `pid_for_rank` maps each neighbour's *rank* to its actual
Distributed pid. The default `r -> r + 1` assumes the standard Julia
convention (master pid 1, auto-spawned workers pid 2, 3, …). Pass a
custom function when the standard convention does not hold (e.g., a
cluster manager assigns non-contiguous pids).

# Rendezvous protocol

`sendto_async` is fire-and-forget; `gather` only guarantees that every
participant *reached* the barrier, not that the asynchronous sendto
tasks scheduled by the peer have *run*. To make the protocol race-free
without polling, the barrier is split into two phases:

1. **Phase 1** — schedule one `sendto_async` per neighbour and
   `fetch` the futures.  Fetching on master blocks master until its
   own sendto task (running on the worker) completes; while master is
   blocked, master's local task queue can process the worker's
   sendto task (which writes to *master's* HaloStash) so the entry
   becomes visible *before* master proceeds to Phase 2.
2. **Phase 2** — `gather` as a "we both finished Phase 1" barrier.
   Both sides have by now populated each other's HaloStashes.
3. **Phase 3** — read `HaloStash.incoming` and apply to the local
   ghost zone.
"""
function exchange_boundary_data(
    mesh::DistributedMesh,
    field_name::Symbol = :rho;
    pid_for_rank::Function = r -> r + 1,
)
    isempty(mesh.neighbors) && return nothing

    # ---- Phase 1: schedule sends, then block on each so the receiver's
    #                task is forced to run while we are waiting.
    send_futures = Future[]
    peer_pids    = Int[]
    for neighbor in mesh.neighbors
        direction = _direction_to(mesh.rank, neighbor, mesh.size, mesh.config.dim)
        halo      = _extract_halo(mesh, direction, field_name)
        peer_pid  = pid_for_rank(neighbor)
        push!(peer_pids, peer_pid)
        # sendto_async(p, expr, data, mod) schedules
        # `Core.eval(mod, Expr(:(=), expr, data))` on worker `p`, so
        # the assignment lands in the *receiver's* HaloStash (which
        # is process-global on that worker).
        f = sendto_async(
            peer_pid,                                              # p
            :(incoming[($(mesh.rank), $(QuoteNode(field_name)))]), # expr (LHS)
            halo,                                                   # data (RHS)
            HaloStash,                                              # mod
        )
        push!(send_futures, f)
    end
    foreach(fetch, send_futures)

    # ---- Phase 2: barrier.  Both sides have finished Phase 1, so the
    #                peer's sendto task that writes to our HaloStash
    #                has already been processed by our local queue.
    all_pids = vcat(peer_pids, myid())
    gather(all_pids, :myid, ParallelOperations)

    # ---- Phase 3: apply received halos ----------------------------------
    for neighbor in mesh.neighbors
        direction = _direction_to(mesh.rank, neighbor, mesh.size, mesh.config.dim)
        key = (neighbor, field_name)
        haskey(HaloStash.incoming, key) || error(
            "HaloStash missing entry for (neighbor=$neighbor, field=$field_name) " *
            "after barrier — neighbour may not have written its slice."
        )
        _apply_halo!(mesh, HaloStash.incoming[key], direction, field_name)
        delete!(HaloStash.incoming, key)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Synchronise every standard field (`rho`, `phi`, `e`) across the
neighbour processes of `mesh`. Add more field names here if your
physics uses them.
"""
function sync_mesh(mesh::DistributedMesh)
    for field_name in (:rho, :phi, :e)
        exchange_boundary_data(mesh, field_name)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Synchronise a single named field across all neighbour processes.
"""
function sync_field(mesh::DistributedMesh, field_name::Symbol)
    return exchange_boundary_data(mesh, field_name)
end
