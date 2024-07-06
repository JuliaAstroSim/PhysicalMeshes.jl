"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Plane{T} <: AbstractPlane{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
end

normal(a,b,c) = normalize(cross(b-a, ustrip(c-a)))
normal(plane::Plane) = normal(plane.a, plane.b, plane.c)

distance(p::AbstractPoint3D, plane::Plane) = abs(dot(normal(plane), p-plane.a))

coplanar(p::AbstractPoint3D, plane::Plane, threshold::Number) = distance(p, plane) <= threshold