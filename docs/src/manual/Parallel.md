```@raw html
<!--
  Parallel backend / operations manual page.
-->
```

# Parallel Backend & Operations

`PhysicalMeshes.jl` carries a small but complete adaptive-parallel layer
that hides whether the underlying storage is a host `Array`, a GPU
`CuArray`, or a distributed `DArray`. The rest of the
JuliaAstroSim ecosystem — `WaveDM.jl`, `PhysicalSPH.jl`,
`PhysicalPlasma.jl` — is built on top of it.

## Backend

See [Library reference](@ref Methods) for the canonical docstrings of
`ParallelBackend`, `detect_resources`, `select_backend`, and the rest of
this page's symbols (`is_serial`, `to_device`, `parallel_sum`, …).

`ParallelBackend` is the value type that selects how to dispatch
work. Its `kind` field is one of `:serial`, `:threads`, `:distributed`,
or `:gpu`; the closures stored on it (`device_array`, `distribute_array`,
`local_array`, `release!`) are picked at construction time.

`detect_resources()` returns a `NamedTuple` of available counts:

```julia
(nthreads = 1, nworkers = 0, pids = Int[], has_cuda = false,
 gpu_count = 0, gpu_name = "", has_darrays = false,
 has_parallelops = false)
```

`select_backend(; gpu = false, distributed = false, kind = :auto)`
chooses the most appropriate backend for the host. The selection rules:

1. `:gpu` if `gpu = true` and CUDA is functional
2. `:distributed` if `distributed = true` and ≥ 2 workers + DistributedArrays
3. `:threads` if more than 1 Julia thread
4. `:serial` otherwise

Passing `kind = :threads` (or any explicit kind) forces the selector to
honour the request; if the requested backend cannot be satisfied, a
warning is emitted and the next-best option is picked.

## Convenience predicates

See [Library reference](@ref Methods) for `is_serial`, `is_threads`,
`is_distributed`, `is_gpu`, `is_parallel`.

## Conversions

See [Library reference](@ref Methods) for `to_device`, `to_host`,
`distribute`, `release!`.

- `to_device(A, backend)` — `cu(A)` for `:gpu`; identity elsewhere.
- `to_host(A, backend)` — `Array(A)` for `:gpu`, `collect(A)` for
  `:distributed`; identity elsewhere.
- `distribute(A, backend)` — `DistributedArrays.distribute(A)` for
  `:distributed`; identity elsewhere.
- `release!(A, backend)` — `CUDA.unsafe_free!(A)` for `:gpu`; no-op
  elsewhere.

## Reductions and parallel primitives

See [Library reference](@ref Methods) for `parallel_sum`, `parallel_maximum`,
`parallel_minimum`, `parallel_findmax`, `parallel_quantile`,
`parallel_sumprod`, `parallel_broadcast!`.

Most of these wrap the corresponding `Base` methods. The notable
exceptions:

- `parallel_findmax` collapses `DArray.findmax`'s two-index tuple
  output to a single global maximum.
- `parallel_quantile(values, weights, qs, backend)` calls
  `Statistics.quantile` with `StatsBase.weights(...)` (see the
  `Statistics` vs `StatsBase` note in the changelog). For `:distributed`
  it gathers first because `quantile` is not defined for `DArray`.

## FFTW thread pinning

See [Library reference](@ref Methods) for `set_backend_fft_threads!`.

`set_backend_fft_threads!(backend)` pins the local-process FFTW thread
count to `backend.fft_threads`. The default per-kind thread counts are:

| `backend.kind`   | `fft_threads`           |
|------------------|-------------------------|
| `:serial`        | 1                       |
| `:threads`       | `Threads.nthreads()`    |
| `:distributed`   | 1 per worker            |
| `:gpu`           | (no-op, FFTW is GPU-side)|

The distributed-backend choice of `1` per worker trades off per-worker
peak FFT speed for predictability — keeping the working set predictable
across workers avoids OOMs when running with full-array replication.

## Why not Dagger.jl?

Dagger.jl's `DRef` finalizer chain does not reliably release remote
memory, which has shown up as monotonically growing DArray allocations
in our benchmarks. The base `DistributedArrays.distribute` (without
Dagger) gives us raw `DArray`s whose lifecycle is fully controlled by
Julia's normal GC, so the parallel backend uses that.

## Example: a minimal backend-aware loop

```julia
using PhysicalMeshes

backend = select_backend()
A = rand(Float64, 1_000_000)
B = similar(A)

parallel_sum(A, backend)             # ≈ 5e5 ± float
parallel_findmax(A, backend)
parallel_broadcast!(x -> x + 1.0, B, A, backend)
```

If `CUDA.functional()` is true and `gpu = true` is passed to
`select_backend`, the same code runs against a `CuArray` without further
changes.
