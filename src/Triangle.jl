struct Triangle2D{T<:Union{Number, Quantity}} <: AbstractTriangle2D{T}
    a::PVector2D{T}
    b::PVector2D{T}
    c::PVector2D{T}
end

struct Triangle3D{T<:Union{Number, Quantity}} <: AbstractTriangle3D{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
end

len(t::AbstractTriangle) = norm(t.a - t.b) + norm(t.a - t.c) + norm(t.b - t.c)