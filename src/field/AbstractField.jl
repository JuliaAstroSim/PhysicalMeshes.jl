abstract type AbstractField{T} end

@inline iterate(f::T) where T <: AbstractField = (f, nothing)
@inline iterate(f::T, st) where T <: AbstractField = nothing

function Base.show(io::IO, f::AbstractField)
    print(io, typeof(f))
end

abstract type ScalarField{T} <: AbstractField{T} end
abstract type VectorField{T} <: AbstractField{T} end
abstract type TensorField{T} <: AbstractField{T} end

@inline getindex(f::AbstractField, i...) = error("getindex not implemented for ", typeof(f))
@inline setindex!(f::AbstractField, v, i...) = error("setindex! not implemented for ", typeof(f))
@inline size(f::AbstractField) = error("size not implemented for ", typeof(f))
@inline length(f::AbstractField) = error("length not implemented for ", typeof(f))
