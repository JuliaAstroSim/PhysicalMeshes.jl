"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Ray2D{T<:Number} <: AbstractRay2D{T}
    "origin of ray"
    x::PVector2D{T}
    "direction vector"
    n::PVector2D{T}
    "azimuth angle"
    theta
end

"The most often used constructor"
Ray2D(x::PVector2D{T}, n::PVector2D{T}) where T<:Number = Ray2D(x, n, atan(n.y, n.x))
Ray2D(n::PVector2D{T}) where T<:Number = Ray2D(PVector2D(T), n, atan(n.y, n.x))
"The most often used constructor"
Ray2D(x::PVector2D{T}, theta::Number) where T<:Number = Ray2D(x, PVector2D(T(sin(theta)), T(cos(theta))), theta)
Ray2D(theta::Number) = Ray2D(PVector2D(), PVector2D(sin(theta), cos(theta)), theta)

zero(r::Ray2D) = Ray2D(r.x, zero(r.n), zero(r.x.x))
iszero(r::Ray2D) = iszero(r.n)

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Ray3D{T<:Number} <: AbstractRay3D{T}
    "origin of ray"
    x::PVector{T}
    "direction vector"
    n::PVector{T}
end

zero(r::Ray3D) = Ray3D(r.x, zero(r.n), zero(r.x.x))
iszero(r::Ray3D) = iszero(r.n)

### Ray2D reflect from Line2D
"""
$(TYPEDSIGNATURES)
"""
function intersect(ray::Ray2D, line::Line2D)
    AB = line.b - line.a
    ray_norm = PVector2D(-ray.n.y, ray.n.x) # ⟂ ray
    dot_ray_norm_AB = dot(ray_norm, AB)
    if iszero(dot_ray_norm_AB)
        return false, PVector2D()  # ∥ ray, no intersection
    end
    AX = ray.x - line.a
    t1 = dot(PVector2D(-AB.y, AB.x), AX) / dot_ray_norm_AB
    t2 = dot(AX, ray_norm) / dot_ray_norm_AB
    if t1 >= 0 && t2 >= 0 && t2 <= one(t2)
        intersection_point = ray.x + ray.n * t1 # rescaling
        return true, intersection_point
    else
        return false, PVector2D()  # no intersection
    end
end

"""
$(TYPEDSIGNATURES)
"""
reflect(vec::AbstractPoint, n::AbstractPoint) = vec - 2*ustrip(n)*dot(vec, ustrip(n))

"""
$(TYPEDSIGNATURES)
"""
function reflect(ray::Ray2D, line::Line2D)
    hit, intersection = intersect(ray, line)
    if hit
        return Ray2D(intersection, reflect(ray.n, normal(line)))
    else
        return nothing
    end
end

### Ray3D reflect from Plane
"""
$(TYPEDSIGNATURES)
"""
function intersect(ray::Ray3D, plane::Plane)
    N = normal(plane)
    denom = dot(ray.n, N)
    if norm(denom) < 1e-6 * unit(denom)
        return false, PVector(0, 0, 0)  # ∥, no intersection
    end
    t = dot(plane.a - ray.x, N) / denom
    if t < 0
        return false, PVector(0, 0, 0)  # no intersection
    end
    intersection_point = ray.x + ray.n * t # rescaling
    return true, intersection_point
end

"""
$(TYPEDSIGNATURES)
"""
function reflect(ray::Ray3D, plane::Plane)
    hit, intersection = intersect(ray, plane)
    if hit
        return Ray3D(intersection, reflect(ray.n, normal(plane)))
    else
        return nothing
    end
end