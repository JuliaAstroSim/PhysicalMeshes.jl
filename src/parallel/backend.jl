# Parallel backend abstraction for the JuliaAstroSim ecosystem.
#
# This module provides an adaptive parallel backend that can transparently
# route work to the most appropriate execution engine:
#
#   :serial       - single thread, no distribution
#   :threads      - multi-threaded BLAS / FFTW / Strided
#   :distributed  - DistributedArrays.jl (DArray) + ParallelOperations.jl
#                   for cross-worker communication
#   :gpu          - CUDA.jl (CuArray)
#
# The backend is selected automatically from the host environment but can
# also be forced.  All numerical kernels (Poisson, FFT, reductions) are
# dispatched on the backend kind so the rest of the ecosystem can stay
# backend-agnostic.
#
# Dagger.jl is intentionally NOT used here: the DRef finalizer chain
# (issue JuliaParallel/Dagger.jl#649) cannot reliably release remote
# memory, and our previous testing confirmed that DArray allocations
# accumulate monotonically when Dagger.jl is in the call path.  We rely
# instead on `DistributedArrays.distribute` (without Dagger) which gives
# us raw DArrays whose lifecycle is fully controlled by Julia's GC.


"""
    ParallelBackend

A value type describing the parallel execution context used by the
JuliaAstroSim packages (`WaveDM.jl`, `PhysicalSPH.jl`,
`PhysicalPlasma.jl`, …).

# Fields
- `kind::Symbol`            One of `:serial`, `:threads`, `:distributed`, `:gpu`
- `gpu::Bool`               Whether GPU acceleration is active
- `nthreads::Int`           Number of Julia threads (1 means serial)
- `nworkers::Int`           Number of Distributed worker processes
- `pids::Vector{Int}`       Worker pids (empty for non-distributed backends)
- `has_darrays::Bool`       Whether DistributedArrays.jl is loaded
- `has_parallelops::Bool`   Whether ParallelOperations.jl is loaded
- `device_array`            Function `A -> cu(A)` (GPU) / `A -> A` (others)
- `distribute_array`        Function `A -> DArray(A, pids, dist)` /
                            `A -> A`
- `local_array`             Function `A -> Array(A)` (collect DArray / CuArray)
- `release!`                In-place device-memory release
- `fft_threads`             Number of FFTW threads (1 for distributed to keep
                            the worker-local FFT cheap and predictable)

Use [`select_backend`](@ref) to construct one.
"""
struct ParallelBackend
    kind::Symbol
    gpu::Bool
    nthreads::Int
    nworkers::Int
    pids::Vector{Int}
    has_darrays::Bool
    has_parallelops::Bool
    device_array        # Function: move host array to device (cu / identity)
    distribute_array    # Function: distribute to workers (DArray / identity)
    local_array         # Function: collect to host (Array / identity)
    release!            # In-place: release device memory (CUDA.unsafe_free!)
    fft_threads::Int
end

function Base.show(io::IO, b::ParallelBackend)
    print(io, "ParallelBackend(", b.kind,
        ", gpu=", b.gpu,
        ", nthreads=", b.nthreads,
        ", nworkers=", b.nworkers, ")")
end

# ----------------------------------------------------------------------------
# Resource detection
# ----------------------------------------------------------------------------

"""
    detect_resources()

Probe the Julia process for available parallel hardware.

Returns a NamedTuple with the available counts.  Callers can use this to
make an informed decision before constructing a backend.
"""
function detect_resources()
    nthreads = Threads.nthreads()
    nworkers = try
        Distributed.nworkers()
    catch
        0
    end
    pids = try
        Distributed.workers()
    catch
        Int[]
    end
    has_cuda = false
    gpu_count = 0
    gpu_name = ""
    if @isdefined(CUDA)
        try
            if CUDA.functional()
                has_cuda = true
                gpu_count = length(CUDA.devices())
                gpu_name = string(CUDA.name(CUDA.device()))
            end
        catch
            # CUDA not initialized or no driver
        end
    end
    has_darrays = (@isdefined(DistributedArrays))
    has_parallelops = (@isdefined(ParallelOperations))
    return (
        nthreads = nthreads,
        nworkers = nworkers,
        pids = pids,
        has_cuda = has_cuda,
        gpu_count = gpu_count,
        gpu_name = gpu_name,
        has_darrays = has_darrays,
        has_parallelops = has_parallelops,
    )
end

# ----------------------------------------------------------------------------
# Backend selection
# ----------------------------------------------------------------------------

"""
    select_backend(; gpu::Bool=false, distributed::Bool=false,
                   kind::Symbol=:auto) -> ParallelBackend

Pick the most appropriate parallel backend for the current host.

# Keyword arguments
- `gpu::Bool`         Hint that GPU should be preferred when available.
- `distributed::Bool` Hint that the simulation should use multiple
                       worker processes when available.
- `kind::Symbol`      `:auto` (default) → pick the best; or one of
                       `:serial`, `:threads`, `:distributed`, `:gpu`.

# Selection rules (when `kind=:auto`)
1. If `gpu=true` and CUDA is functional → `:gpu`
2. If `distributed=true` and at least 2 workers → `:distributed`
3. If more than 1 Julia thread → `:threads`
4. Otherwise → `:serial`

When a requested backend cannot be satisfied (e.g. `:gpu` requested but
no CUDA), the selector falls back to the next-best option and logs a
warning.
"""
function select_backend(; gpu::Bool=false, distributed::Bool=false, kind::Symbol=:auto)
    res = detect_resources()

    chosen = kind
    if chosen === :auto
        if gpu && res.has_cuda
            chosen = :gpu
        elseif distributed && res.nworkers >= 2 && res.has_darrays
            chosen = :distributed
        elseif res.nthreads > 1
            chosen = :threads
        else
            chosen = :serial
        end
    end

    # Validate the chosen kind against the environment, fall back if needed.
    if chosen === :gpu && !res.has_cuda
        @warn "GPU backend requested but CUDA is not functional; falling back to threads"
        chosen = res.nthreads > 1 ? :threads : :serial
    elseif chosen === :distributed && (res.nworkers < 2 || !res.has_darrays)
        @warn "Distributed backend requested but workers<2 or DistributedArrays not loaded; falling back to threads"
        chosen = res.nthreads > 1 ? :threads : :serial
    elseif chosen === :threads && res.nthreads <= 1
        @warn "Threads backend requested but only 1 thread; falling back to serial"
        chosen = :serial
    end
    chosen in (:serial, :threads, :distributed, :gpu) ||
        throw(ArgumentError("Unknown backend kind: $(chosen)"))

    return _build_backend(chosen, res)
end

# ----------------------------------------------------------------------------
# Backend construction
# ----------------------------------------------------------------------------

function _build_backend(kind::Symbol, res)
    gpu = (kind === :gpu)
    # Pick FFT threads: distributed backend uses 1 thread per worker to keep
    # memory pressure predictable (the whole array lives on each worker);
    # threads/gpu use all available threads; serial uses 1.
    fft_threads = if kind === :distributed
        1
    elseif kind === :threads
        res.nthreads
    else
        1
    end

    device_array = _make_device_array(kind, res)
    distribute_array = _make_distribute_array(kind, res)
    local_array = _make_local_array(kind)
    release! = _make_release(kind)
    pids = kind === :distributed ? collect(res.pids) : Int[]

    return ParallelBackend(
        kind, gpu,
        res.nthreads, res.nworkers, pids,
        res.has_darrays, res.has_parallelops,
        device_array, distribute_array, local_array, release!,
        fft_threads,
    )
end

# ---- helpers: per-kind function factories ---------------------------------

# Move a host array to the device memory of this backend
function _make_device_array(kind::Symbol, res)
    if kind === :gpu
        return A -> CUDA.cu(A)
    else
        # For non-GPU, "device" = host.  The identity function preserves
        # the array as-is (Array, DArray, SubArray, ...).
        return identity
    end
end

# Distribute an array across workers
function _make_distribute_array(kind::Symbol, res)
    if kind === :distributed
        nw = max(1, res.nworkers)
        pids = collect(res.pids)
        return A -> DistributedArrays.distribute(A; procs=pids, dist=[1, 1, nw])
    else
        return identity
    end
end

# Collect a (potentially distributed / device-resident) array to host
function _make_local_array(kind::Symbol)
    if kind === :distributed
        return A -> collect(A)
    elseif kind === :gpu
        return A -> Array(A)
    else
        return identity
    end
end

# Release device memory (no-op for serial / distributed)
function _make_release(kind::Symbol)
    if kind === :gpu
        return A -> (CUDA.unsafe_free!(A); nothing)
    else
        return A -> nothing
    end
end

# ----------------------------------------------------------------------------
# Convenience predicates
# ----------------------------------------------------------------------------

is_serial(b::ParallelBackend)       = b.kind === :serial
is_threads(b::ParallelBackend)      = b.kind === :threads
is_distributed(b::ParallelBackend)  = b.kind === :distributed
is_gpu(b::ParallelBackend)          = b.kind === :gpu
is_parallel(b::ParallelBackend)     = b.kind !== :serial
