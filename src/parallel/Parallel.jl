include("ParallelBackend.jl")
include("DistributedMesh.jl")
include("Partitioning.jl")
include("Communication.jl")

export
    # Parallel backend (adaptive: :serial / :threads / :distributed / :gpu)
    ParallelBackend,
    select_backend, detect_resources,
    is_serial, is_threads, is_distributed, is_gpu, is_parallel,
    set_backend_fft_threads!,
    to_device, to_host, distribute, release!,
    parallel_sum, parallel_maximum, parallel_minimum,
    parallel_findmax, parallel_quantile, parallel_sumprod,
    parallel_broadcast!,

    # Mesh / field distribution
    DistributedMesh,
    PartitioningStrategy,
    UniformPartitioning,
    BalancedPartitioning,
    partition_mesh,
    exchange_boundary_data,
    sync_mesh,
    sync_field
