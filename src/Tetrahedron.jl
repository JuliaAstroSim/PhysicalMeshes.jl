struct Tetrahedron{T<:Number} <: AbstractTetrahedron{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
    d::PVector{T}
end

