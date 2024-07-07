"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Polygon2D{T} <: AbstractPolygon2D{T}
    vertices::AbstractVector{T}

    function Polygon2D(vertices)
        if length(vertices) < 3
            error("Number of vertices less than 3, a concrete polygon should have more than 2 vertices")
        end
    
        new{eltype(vertices)}(vertices)
    end
end

"""
$(TYPEDSIGNATURES)
Compute the unit normal vector from the first three vertices
"""
function normal(polygon::Polygon2D)
    x = typeof(polygon.vertices[1].x)
    PVector(zero(x), zero(x), one(x) * unit(x)) # In Unitful.jl, `zero` preserves units whereas `one` does not.
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Polygon3D{T} <: AbstractPolygon3D{T}
    vertices::AbstractVector{T}

    function Polygon3D(vertices)
        if length(vertices) < 3
            error("Number of vertices less than 3, a concrete polygon should have more than 2 vertices")
        end
    
        if !coplanar(vertices, unit(vertices[1].x)^2 * 1e-6)
            @warn "The polygon is not coplanar with 1e-6 threshold!"
        end
    
        new{eltype(vertices)}(vertices)
    end
end

"""
$(TYPEDSIGNATURES)
Compute the unit normal vector from the first three vertices
"""
normal(polygon::Polygon3D) = normal(polygon.vertices[1:3]...)

"""
$(TYPEDSIGNATURES)
The internal angles of a convex polygon are all smaller than π,
or equivalently, the cross products of any two adjacent edges have the same direction.
"""
function isconvex(polygon::AbstractPolygon)
    n = length(polygon.vertices)
    if n < 4
        return true
    end

    N = normal(polygon)
    cross_sign = 0
    for i in 1:n
        v1 = polygon.vertices[i]
        v2 = polygon.vertices[mod1(i+1, n)]
        v3 = polygon.vertices[mod1(i+2, n)]
        if cross_sign == 0
            cross_sign = sign(dot(cross(v2 - v1, v3 - v2), N))
        else
            if sign(dot(cross(v2 - v1, v3 - v2), N)) != cross_sign
                return false
            end
        end
    end
    return true
end

"""
$(TYPEDSIGNATURES)

First construct a plane from the first three vertices, then iteratively check whether all other vertices are coplanar with this plane or not.
"""
function coplanar(vertices, threshold)
    n = length(vertices)
    
    # A concrete polygon have more than 2 vertices
    if n == 3
        return true
    else
        plane = Plane(vertices[1:3]...)
        for i = 4:n
            if !coplanar(vertices[i], plane, threshold)
                return false
            end
        end
        return true
    end
end

coplanar(polygon::Polygon3D, threshold) = coplanar(polygon.vertices, threshold)

"""
$(TYPEDSIGNATURES)
For convex polygons, use ray casting method to check whether the point is inside the polygon.
"""
function is_inbound_ray_casting(pos::PVector2D, polygon::Polygon2D)
    vertices = polygon.vertices
    n = length(vertices)
    inside = false
    j = n
    for i in 1:n
        if (vertices[i].y > pos.y) != (vertices[j].y > pos.y) &&
           (pos.x < (vertices[j].x - vertices[i].x) * (pos.y - vertices[i].y) / (vertices[j].y - vertices[i].y) + vertices[i].x)
            inside = !inside
        end
        j = i
    end
    return inside
end

#TODO ray-casting method for Polygon3D

"""
$(TYPEDSIGNATURES)
For convex polygons, use cross product method to check whether the point is inside the polygon.
"""
function is_inbound_cross_product(pos, polygon)
    vertices = polygon.vertices
    n = length(vertices)
    N = normal(polygon)
    cross_sign = 0
    for i in 1:n
        v1 = vertices[i]
        v2 = vertices[mod1(i+1, n)]
        edge = v2 - v1
        vp = pos - v1
        if cross_sign == 0
            cross_sign = sign(dot(cross(edge, vp), N))
        else
            if sign(dot(cross(edge, vp), N)) != cross_sign
                return false
            end
        end
    end
    return true
end

function is_inbound(pos, polygon)
    if isconvex(polygon)
        return is_inbound_cross_product(pos, polygon)
    else
        return is_inbound_ray_casting(pos, polygon)
    end
end

### Some common polygons

function polygon_rect()
    
end

function polygon_regular()
    
end