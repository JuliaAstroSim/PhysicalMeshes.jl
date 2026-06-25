struct ArrayScalarField{T, A <: AbstractArray{T}} <: ScalarField{T}
    data::A
end

struct ArrayVectorField{T, A <: AbstractArray{T}} <: VectorField{T}
    data::A
    dim::Int
end

struct ArrayTensorField{T, A <: AbstractArray{T}} <: TensorField{T}
    data::A
    dim::Int
end

@inline getindex(f::ArrayScalarField, i...) = f.data[i...]
@inline setindex!(f::ArrayScalarField, v, i...) = f.data[i...] = v
@inline size(f::ArrayScalarField) = size(f.data)
@inline size(f::ArrayScalarField, d::Int) = size(f.data, d)
@inline length(f::ArrayScalarField) = length(f.data)
@inline ndims(f::ArrayScalarField) = ndims(f.data)
@inline ndims(::Type{ArrayScalarField{T,A}}) where {T,A} = ndims(A)
@inline axes(f::ArrayScalarField) = axes(f.data)
@inline axes(f::ArrayScalarField, d::Int) = axes(f.data, d)
@inline eltype(f::ArrayScalarField{T}) where T = T
@inline eltype(::Type{ArrayScalarField{T,A}}) where {T,A} = T
@inline Base.lastindex(f::ArrayScalarField) = lastindex(f.data)
@inline Base.lastindex(f::ArrayScalarField, d::Int) = lastindex(f.data, d)
@inline Base.firstindex(f::ArrayScalarField) = firstindex(f.data)
@inline Base.firstindex(f::ArrayScalarField, d::Int) = firstindex(f.data, d)

@inline function getindex(f::ArrayVectorField, i...)
    n_spatial = ndims(f.data) - 1
    if length(i) == ndims(f.data)
        # Caller indexed into the underlying data directly, including the
        # trailing vector-component dim. This happens in 3D FDTD loops
        # like `B[i, j, k, 3]` where the user actually wants the raw
        # scalar at that index, not the vector slice.
        return f.data[i...]
    end
    idx = CartesianIndices(size(f.data)[1:n_spatial])[i...]
    return f.data[idx, :]
end

@inline function setindex!(f::ArrayVectorField, v, i...)
    n_spatial = ndims(f.data) - 1
    if length(i) == ndims(f.data)
        return setindex!(f.data, v, i...)
    end
    idx = CartesianIndices(size(f.data)[1:n_spatial])[i...]
    f.data[idx, :] .= v
end

@inline size(f::ArrayVectorField) = size(f.data)[1:end-1]
@inline size(f::ArrayVectorField, d::Int) = size(f.data, d)
@inline length(f::ArrayVectorField) = prod(size(f))
@inline ndims(f::ArrayVectorField) = ndims(f.data) - 1
@inline ndims(::Type{ArrayVectorField{T,A}}) where {T,A} = ndims(A) - 1
@inline axes(f::ArrayVectorField) = axes(f.data)[1:end-1]
@inline axes(f::ArrayVectorField, d::Int) = axes(f.data, d)
@inline eltype(f::ArrayVectorField{T}) where T = T
@inline eltype(::Type{ArrayVectorField{T,A}}) where {T,A} = T
# lastindex/firstindex for ArrayVectorField: spatial dims only (vector dim stripped)
# `f[end]` must resolve to the last *spatial* element; `f.data` has an extra
# trailing dim holding the vector components.
@inline Base.lastindex(f::ArrayVectorField) = length(f)
@inline Base.lastindex(f::ArrayVectorField, d::Int) = size(f.data, d)
@inline Base.firstindex(f::ArrayVectorField) = 1
@inline Base.firstindex(f::ArrayVectorField, d::Int) = 1

@inline function getindex(f::ArrayTensorField, i...)
    idx = CartesianIndices(size(f.data)[1:end-2])[i...]
    return reshape(f.data[idx, :, :], (f.dim, f.dim))
end

@inline function setindex!(f::ArrayTensorField, v, i...)
    idx = CartesianIndices(size(f.data)[1:end-2])[i...]
    f.data[idx, :, :] .= v
end

@inline size(f::ArrayTensorField) = size(f.data)[1:end-2]
@inline size(f::ArrayTensorField, d::Int) = size(f.data, d)
@inline length(f::ArrayTensorField) = prod(size(f))
@inline ndims(f::ArrayTensorField) = ndims(f.data) - 2
@inline ndims(::Type{ArrayTensorField{T,A}}) where {T,A} = ndims(A) - 2
@inline axes(f::ArrayTensorField) = axes(f.data)[1:end-2]
@inline axes(f::ArrayTensorField, d::Int) = axes(f.data, d)
@inline eltype(f::ArrayTensorField{T}) where T = T
@inline eltype(::Type{ArrayTensorField{T,A}}) where {T,A} = T
# lastindex/firstindex for ArrayTensorField: spatial dims only (last two trailing
# dims hold the (dim x dim) tensor components and are not exposed by indexing).
@inline Base.lastindex(f::ArrayTensorField) = length(f)
@inline Base.lastindex(f::ArrayTensorField, d::Int) = size(f.data, d)
@inline Base.firstindex(f::ArrayTensorField) = 1
@inline Base.firstindex(f::ArrayTensorField, d::Int) = 1

function ArrayScalarField(::Type{T}, dims::NTuple{N, Int}) where {T, N}
    data = zeros(T, dims)
    return ArrayScalarField(data)
end

function ArrayVectorField(::Type{T}, dims::NTuple{N, Int}, dim::Int) where {T, N}
    data = zeros(T, dims..., dim)
    return ArrayVectorField(data, dim)
end

function ArrayTensorField(::Type{T}, dims::NTuple{N, Int}, dim::Int) where {T, N}
    data = zeros(T, dims..., dim, dim)
    return ArrayTensorField(data, dim)
end

# Numerical reductions and arithmetic on ArrayScalarField so that
# expressions like `sum(m.rho) * Δ^3` and `m.rho .* m.pos` work
# without having to reach for `m.rho.data` explicitly.
@inline Base.sum(f::ArrayScalarField; kw...) = sum(f.data; kw...)
@inline Base.sum(f::ArrayVectorField; kw...) = sum(f.data; kw...)
@inline Base.prod(f::ArrayScalarField; kw...) = prod(f.data; kw...)
@inline Base.minimum(f::ArrayScalarField; kw...) = minimum(f.data; kw...)
@inline Base.maximum(f::ArrayScalarField; kw...) = maximum(f.data; kw...)
@inline Base.extrema(f::ArrayScalarField; kw...) = extrema(f.data; kw...)

@inline Base.:*(f::ArrayScalarField, x::Number) = ArrayScalarField(f.data .* x)
@inline Base.:*(x::Number, f::ArrayScalarField) = ArrayScalarField(f.data .* x)
@inline Base.:/(f::ArrayScalarField, x::Number) = ArrayScalarField(f.data ./ x)
@inline Base.:+(f::ArrayScalarField, x::Number) = ArrayScalarField(f.data .+ x)
@inline Base.:+(x::Number, f::ArrayScalarField) = ArrayScalarField(f.data .+ x)
@inline Base.:-(f::ArrayScalarField, x::Number) = ArrayScalarField(f.data .- x)

# Element-wise product between a scalar field and an array of positions
# (StructArray of PVector) so `m.rho .* m.pos` works.
@inline function Base.broadcasted(*, f::ArrayScalarField, sa::StructArray)
    return f.data .* sa
end

# Broadcast assignment support: `m.rho .= expr` calls materialize! -> copyto!
@inline function Base.copyto!(f::ArrayScalarField, bc::Base.Broadcast.Broadcasted)
    copyto!(f.data, bc)
    return f
end
@inline function Base.copyto!(f::ArrayVectorField, bc::Base.Broadcast.Broadcasted)
    copyto!(f.data, bc)
    return f
end

# Broadcast destination for slice assignment `f[I...] .= value`:
# `getindex` returns a (possibly copied) value, so without `dotview` the
# broadcast would write into the copy and the slice assignment would be
# a silent no-op on `f.data`. Forward to a real `view` of `f.data` so
# `m.phi[2:end-1, ...] .= rhs` (the pattern used in PhysicalFFT.jl:342)
# actually updates `m.phi.data`.
@inline Base.Broadcast.dotview(f::ArrayScalarField, i...) = Base.view(f.data, i...)

# broadcastable: let broadcast use the underlying data array
@inline Base.broadcastable(f::ArrayScalarField) = f.data
@inline Base.broadcastable(f::ArrayVectorField) = f.data
@inline Base.broadcastable(f::ArrayTensorField) = f.data

function Base.show(io::IO, f::ArrayScalarField)
    print(io, "ArrayScalarField{$(eltype(f.data))} with size $(size(f.data))")
end

function Base.show(io::IO, f::ArrayVectorField)
    print(io, "ArrayVectorField{$(eltype(f.data))} with size $(size(f.data)[1:end-1]) and dim $(f.dim)")
end

function Base.show(io::IO, f::ArrayTensorField)
    print(io, "ArrayTensorField{$(eltype(f.data))} with size $(size(f.data)[1:end-2]) and dim $(f.dim)x$(f.dim)")
end
