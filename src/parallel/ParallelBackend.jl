# ParallelBackend: adaptive parallel execution for the JuliaAstroSim
# ecosystem (WaveDM.jl, PhysicalSPH.jl, PhysicalPlasma.jl, …).
#
# Note: this file deliberately lives in `PhysicalMeshes.jl` because the
# mesh / field / partition code is the most natural host for an adaptive
# parallel abstraction (mesh partitioning, boundary exchange, field
# distribution all flow through here).  Downstream packages
# (`WaveDM.jl`, `PhysicalSPH.jl`, `PhysicalPlasma.jl`) `using` this
# module to obtain the same backend-aware primitives.

include("backend.jl")
include("operations.jl")

export
    ParallelBackend,
    select_backend, detect_resources,
    is_serial, is_threads, is_distributed, is_gpu, is_parallel,
    set_backend_fft_threads!,
    to_device, to_host, distribute, release!,
    parallel_sum, parallel_maximum, parallel_minimum,
    parallel_findmax, parallel_quantile, parallel_sumprod,
    parallel_broadcast!
