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
    if t1 >= 0 && t2 >= 0 && t2 <= 1
        intersection_point = ray.x + ray.n * t1
        return true, intersection_point
    else
        return false, PVector2D()  # no intersection
    end
end

reflect(vec::PVector2D, n::PVector2D) = vec - 2*ustrip(n)*dot(vec, ustrip(n))

function reflect(ray::Ray2D, line::Line2D)
    hit, intersection = intersect(ray, line)
    if hit
        return Ray2D(intersection, reflect(ray.n, normal(line)))
    else
        return nothing
    end
end