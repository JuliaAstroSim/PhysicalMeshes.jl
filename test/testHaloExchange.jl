# test/testHaloExchange.jl
#
# Toy test for `exchange_boundary_data` between two ranks.
# Run from the repository root with:
#
#     julia --project test/testHaloExchange.jl
#
# Topology
# --------
# master (pid 1, rank 0) ── worker (pid 2, rank 1)
#         |                                |
#     fills rho = 1.0                fills rho = 2.0
#         |                                |
#         +────── exchange_boundary_data ──+
#
# After exchange:
#   rank 0's right ghost zone (indices 7:8) should hold 2.0
#   rank 1's left  ghost zone (indices 1:2) should hold 1.0
#
# The interior of each rank is left untouched, and that is also asserted.

using Test, Distributed

# Standard Julia convention: master is pid 1, auto-spawned/added workers
# are pid 2, 3, …
#
# Quirk in some Julia 1.12 setups: master pid 1 appears in `workers()`
# initially (with `nworkers() == 1`), so we cannot gate on
# `nworkers() < 1`.  Instead, count "real" workers by filtering out the
# master pid and use `addprocs(1)` only when none exist.
#
# On Windows, `addprocs(1)` can also fail silently if the worker can't
# connect back to the master (firewall blocking the TCP port is the
# usual suspect).  We probe explicitly and emit a clear error.
real_workers = filter(p -> p != 1, workers())
if length(real_workers) < 1
    added = addprocs(1)
    if length(added) < 1
        error(
            "addprocs(1) returned no workers. " *
            "On Windows this usually means the worker process could " *
            "not connect back to the master (firewall, ssh tunnel, " *
            "missing julia on PATH, etc.)."
        )
    end
    real_workers = filter(p -> p != 1, workers())
end

@assert myid() == 1 "Master must be pid 1, got $(myid())"
@assert length(real_workers) >= 1 "Need a real worker (pid != 1); workers() = $(workers())"
peer_pid = real_workers[1]
println("[halo-exchange] master=$(myid()) workers=$(workers()) peer=$(peer_pid)")

@everywhere using PhysicalMeshes

@testset "Halo exchange: 2 ranks, 1D, scalar rho" begin
    NG          = 2            # ghost thickness on each side
    N_local     = 4            # interior cells per rank
    N_with_ghost = N_local + 2 * NG   # = 8 = Len in CellMode

    # pid_for_rank maps a neighbour's rank to its actual pid.  In the
    # standard setup master is rank 0 (pid 1) and the peer is rank 1
    # (pid = peer_pid, which on Julia 1.12 may not be 2).
    function pid_for_rank(r)
        r == 0 ? 1 : peer_pid
    end

    # Per-process setup + exchange.  Returns a NamedTuple describing
    # the local state *after* the halo exchange, so the assertions
    # can run on the master only.
    @everywhere function _local_run(NG, N_with_ghost, peer_pid)
        rank = myid() == 1 ? 0 : 1      # master → 0, peer → 1
        pfr  = r -> r == 0 ? 1 : peer_pid
        # `MeshConfig.Nx` is the interior cell count, so we reconstruct
        # it from the total length.  We deliberately do NOT close over
        # the outer testset's `N_local` variable, because `@everywhere`
        # cannot capture a closure variable from a remote eval context.
        N_local = N_with_ghost - 2 * NG

        cfg = MeshConfig(;
            mode = CellMode(),
            Nx = N_local, Ny = 1, Nz = 1,
            xMin = 0.0, xMax = 1.0,
            NG  = NG,
            dim = 1,
        )
        local_mesh = MeshCartesianStatic(cfg, nothing, nothing)
        fill_value = rank == 0 ? 1.0 : 2.0
        local_mesh.rho.data .= fill_value

        # Build DistributedMesh directly, bypassing the convenience
        # constructor that computes its own local config.  For a 2-rank
        # 1-D chain the neighbour sets are obvious.
        neighbors = rank == 0 ? Int[1] : Int[0]
        mesh = DistributedMesh(cfg, local_mesh, rank, 2, neighbors)

        interior_before = collect(local_mesh.rho.data[NG+1 : N_with_ghost - NG])

        exchange_boundary_data(mesh, :rho; pid_for_rank = pfr)

        data = local_mesh.rho.data
        if rank == 0
            return (
                rank            = rank,
                ghost_side      = :right,
                ghost           = collect(data[N_with_ghost - NG + 1 : N_with_ghost]),
                expected        = 2.0,
                interior_after  = collect(data[NG+1 : N_with_ghost - NG]),
                interior_before = interior_before,
            )
        else
            return (
                rank            = rank,
                ghost_side      = :left,
                ghost           = collect(data[1:NG]),
                expected        = 1.0,
                interior_after  = collect(data[NG+1 : N_with_ghost - NG]),
                interior_before = interior_before,
            )
        end
    end

    # IMPORTANT: spawn worker FIRST.  If master runs _local_run before
    # scheduling the worker, master blocks in Phase 1 (`fetch` of its own
    # sendto_async future) and never reaches the line that would spawn
    # the worker.  Spawning the worker up-front lets it run concurrently
    # with master, so its `sendto_async` task is scheduled onto master's
    # queue while master is blocked, and master's HaloStash is populated
    # before master reaches Phase 3.
    worker_future = @spawnat peer_pid _local_run(NG, N_with_ghost, peer_pid)
    r0 = _local_run(NG, N_with_ghost, peer_pid)                              # master, rank 0
    r1 = fetch(worker_future)                                                # peer, rank 1

    @testset "rank 0 (filled with 1.0) → right ghost = 2.0" begin
        @test r0.rank == 0
        @test r0.ghost_side == :right
        @test length(r0.ghost) == NG
        @test all(r0.ghost .== r0.expected)
        @test r0.interior_after == r0.interior_before   # interior untouched
        @test r0.interior_after == fill(1.0, N_local)
    end

    @testset "rank 1 (filled with 2.0) → left ghost = 1.0" begin
        @test r1.rank == 1
        @test r1.ghost_side == :left
        @test length(r1.ghost) == NG
        @test all(r1.ghost .== r1.expected)
        @test r1.interior_after == r1.interior_before
        @test r1.interior_after == fill(2.0, N_local)
    end

    # Human-friendly echo so the result is visible without `-v`.
    println("\n--- Halo exchange summary ---")
    println("rank 0  interior=$(r0.interior_after)  right ghost=$(r0.ghost)")
    println("rank 1  interior=$(r1.interior_after)  left  ghost=$(r1.ghost)")
end
