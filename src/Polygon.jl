"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Polygon2D{T} <: AbstractPolygon2D{T}
    vertices
end

function Polygon2D(vertices)
    if length(vertices) < 3
        error("Number of vertices less than 3, a concrete polygon should have more than 2 vertices")
    end

    Polygon2D(vertices)
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Polygon3D{T} <: AbstractPolygon3D{T}
    vertices
end

function Polygon3D(vertices)
    if length(vertices) < 3
        error("Number of vertices less than 3, a concrete polygon should have more than 2 vertices")
    end

    if !coplanar(vertices, one(vertices[1].x) * 1e-6)
        @warn "The polygon is not coplanar with 1e-6 threshold!"
    end

    Polygon3D(vertices)
end

function isconvex(polygon::Polygon2D)
    
end

"""
$(TYPEDSIGNATURES)
The internal angles of a convex polygon are all smaller than π,
or equivalently, the cross products of any two adjacent edges have the same direction.
"""
function isconvex(polygon::Polygon3D)
    n = length(polygon.vertices)
    if n < 4
        return true
    end

    N = normal(polygon)
    for i in 1:n
        v1 = polygon.vertices[i]
        v2 = polygon.vertices[mod1(i+1, n)]
        v3 = polygon.vertices[mod1(i+2, n)]
        if dot(cross(v2 - v1, v3 - v2), N) < 0
            return false
        end
    end
    return true
end

"""
$(TYPEDSIGNATURES)
Compute the unit normal vector from the first three vertices
"""
normal(polygon::Polygon3D) = normal(polygon.vertices[1:3]...)

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

function is_inbound(pos, polygon)

end

### Some common polygons

function polygon_rect()
    
end

function polygon_regular()
    
end