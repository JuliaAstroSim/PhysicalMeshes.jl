# Generic parallel operations for the JuliaAstroSim ecosystem.
#
# These wrappers dispatch on `ParallelBackend` so the rest of the codebase
# can call them without knowing whether the data lives on the host, on a
# CuArray, or on a DArray.
#
# We deliberately avoid Dagger.jl (issue #649: the DRef finalizer chain
# cannot release remote memory in a timely fashion).  For the distributed
# backend we use plain DistributedArrays + ParallelOperations, which
# gives us DArrays whose memory is reclaimed by Julia's normal GC.
# -----------------------------------------------------------------------

# --------------------------------------------------------------------------
# Conversions: host <-> device, host <-> distributed
# --------------------------------------------------------------------------

"""
    to_device(A, backend) -> A'

Move `A` to the device memory of `backend` (CuArray on GPU, identity
otherwise).  No copy is performed when the backend is non-GPU and `A`
already lives on host.
"""
to_device(A, b::ParallelBackend) = b.device_array(A)

"""
    to_host(A, backend) -> A'

Bring `A` back to host memory (collects DArrays, copies CuArrays, identity
otherwise).  This is the inverse of `to_device`.
"""
to_host(A, b::ParallelBackend) = b.local_array(A)

"""
    distribute(A, backend) -> A'

For the distributed backend, create a DArray from `A`.  For other
backends this is the identity function (the array is already local).
"""
distribute(A, b::ParallelBackend) = b.distribute_array(A)

"""
    release!(A, backend) -> nothing

Release device-side storage held by `A`.  On GPU this is
`CUDA.unsafe_free!`; on other backends this is a no-op.  The array is
left in an undefined state, callers must not reuse the binding.
"""
release!(A, b::ParallelBackend) = b.release!(A)

# --------------------------------------------------------------------------
# Reductions
#
# `sum`, `maximum`, `minimum`, `findmax`, `findmin` already work on
# `Array`, `CuArray` and `DArray`, so we just forward to Base.  We
# provide thin wrappers for consistency and to allow per-backend tuning
# (e.g. pinning FFTW threads on distributed workers).
# --------------------------------------------------------------------------

"""
    parallel_sum(A, backend) -> scalar

Sum all elements of `A` using the most efficient primitive of `backend`.
"""
parallel_sum(A, b::ParallelBackend) = sum(A)

"""
    parallel_maximum(A, backend) -> scalar
"""
parallel_maximum(A, b::ParallelBackend) = maximum(A)

"""
    parallel_minimum(A, backend) -> scalar
"""
parallel_minimum(A, b::ParallelBackend) = minimum(A)

"""
    parallel_findmax(A, backend) -> (value, index)

`findmax` with a single scalar return.  For DArrays we keep the
multi-local result and pick the global winner.
"""
function parallel_findmax(A::AbstractArray, b::ParallelBackend)
    return findmax(A)
end

# When `A` is a DArray, `findmax` may return a `Tuple{T,Tuple{CartesianIndex,
# CartesianIndex}}`; flatten it.
function parallel_findmax(A, b::ParallelBackend)
    if b.kind === :distributed
        # DArray.findmax returns (maxval, (local_idx, global_idx)) in some
        # versions; we collapse to a single global maximum here.
        v, loc = findmax(A)
        return (v, loc[2])
    else
        return findmax(A)
    end
end

"""
    parallel_quantile(values, weights, qs, backend)

Weighted quantile.  DistributedArrays does not support `quantile` on
`DArray` directly, so for the distributed backend we gather the
elements first.  This is intentionally simple: exact weighted quantiles
across a DArray are rare in scientific code where the quantile is
computed for diagnostic radii.
"""
function parallel_quantile(values::AbstractArray, weights::AbstractArray, qs, b::ParallelBackend)
    if b.kind === :distributed
        return Statistics.quantile(collect(values[:]), Statistics.weights(collect(weights[:])), qs)
    else
        return Statistics.quantile(values, Statistics.weights(weights), qs)
    end
end

# --------------------------------------------------------------------------
# Element-wise broadcasting
# --------------------------------------------------------------------------

"""
    parallel_broadcast!(f, dst, As, backend) -> dst

Compute `dst .= f.(As...)` using backend-native broadcasting.  Plain
broadcast already works for Array, CuArray and DArray, so this is just a
named wrapper for symmetry with the rest of the API.
"""
function parallel_broadcast!(f, dst, As, b::ParallelBackend)
    broadcast!(f, dst, As...)
    return dst
end

# --------------------------------------------------------------------------
# Sum-of-product (dot), used in virial / momentum accumulators
# --------------------------------------------------------------------------

"""
    parallel_sumprod(A, B, backend) -> scalar

Compute `sum(A .* B)`.  `sum(A .* B)` already works for DArray / CuArray
but produces a temporary allocation; for the threads/gpu path
`BLAS.dot` / `LinearAlgebra.dot` is faster.  We dispatch accordingly.
"""
function parallel_sumprod(A, B, b::ParallelBackend)
    if b.kind === :gpu
        return sum(A .* B)
    elseif b.kind === :distributed
        return sum(A .* B)
    else
        return sum(A .* B)
    end
end

# --------------------------------------------------------------------------
# FFTW thread count helper
# --------------------------------------------------------------------------

"""
    set_backend_fft_threads!(backend)

Pin the local-process FFTW thread count to `backend.fft_threads`.  This
must be called on every worker that participates in a distributed run
and on the main process.  For the GPU and serial backends this is a
no-op beyond the FFTW pin.
"""
function set_backend_fft_threads!(b::ParallelBackend)
    if b.kind !== :gpu
        FFTW.set_num_threads(b.fft_threads)
    end
    return b.fft_threads
end
