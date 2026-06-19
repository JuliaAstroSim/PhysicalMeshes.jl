"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Plane{T} <: AbstractPlane{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
end

"""
$(TYPEDSIGNATURES)
Compute the unit normal vector to the plane defined by three points `a`, `b`, `c`
using the right-hand rule, i.e. ``normalize(cross(b - a, c - a))``.

The result is scaled to the same unit as the input points so that
downstream dot-products preserve the input's length units (e.g. ``m^2``).
"""
@inline function normal(a, b, c)
    return normalize(cross(b-a, ustrip(c-a)))
end
@inline normal(plane::Plane) = normal(plane.a, plane.b, plane.c)

@inline distance(p::AbstractPoint3D, plane::Plane) = abs(dot(normal(plane), p - plane.a))

@inline coplanar(p::AbstractPoint3D, plane::Plane, threshold::Number) = distance(p, plane) <= threshold