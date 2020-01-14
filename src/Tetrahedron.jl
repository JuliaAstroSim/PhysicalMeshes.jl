struct Tetrahedron{T<:Union{Number, Quantity}} <: AbstractTetrahedron{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
    d::PVector{T}
end

