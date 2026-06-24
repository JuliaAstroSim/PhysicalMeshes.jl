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

# Compatibility shim: legacy tests and downstream code use
# `polygon.points` while the canonical field is `polygon.vertices`.
@inline Base.getproperty(polygon::Polygon2D, name::Symbol) =
    name === :points ? Base.getfield(polygon, :vertices) : Base.getfield(polygon, name)

@inline Base.getproperty(polygon::Polygon3D, name::Symbol) =
    name === :points ? Base.getfield(polygon, :vertices) : Base.getfield(polygon, name)

@inline Base.propertynames(::Polygon2D) = (:vertices, :points)
@inline Base.propertynames(::Polygon3D) = (:vertices, :points)

"""
$(TYPEDSIGNATURES)
Compute the unit normal vector from the first three vertices
"""
normal(polygon::Polygon3D) = normal(polygon.vertices[1:3]...)

"""
$(TYPEDSIGNATURES)
Compute the unsigned area of a 2D polygon using the shoelace formula.
The polygon is split into triangles anchored at the first vertex for
numerical robustness.
"""
function area(polygon::Polygon2D)
    vs = polygon.vertices
    n = length(vs)
    n < 3 && return zero(promote_type(typeof(vs[1].x), typeof(vs[1].y)))
    acc = zero(promote_type(typeof(vs[1].x), typeof(vs[1].y)))
    @inbounds for i in 2:(n - 1)
        v1 = vs[i] - vs[1]
        v2 = vs[i + 1] - vs[1]
        acc += v1.x * v2.y - v2.x * v1.y
    end
    return abs(acc) / 2
end

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

"""
$(TYPEDSIGNATURES)
Project a 3D point onto the polygon plane and return 2D coordinates in the plane's local coordinate system.
"""
function project_to_plane(pos::PVector, polygon::Polygon3D)
    vertices = polygon.vertices
    N = normal(polygon)
    
    # Create a local coordinate system on the plane
    # Use the first edge as the x-axis. Use ``/norm`` (rather than
    # ``normalize``) so the unit vector is truly dimensionless even
    # for ``Unitful.Quantity`` inputs.
    v0 = vertices[1]
    v1 = vertices[2]
    x_axis = (v1 - v0) / norm(v1 - v0)
    y_axis = cross(N, x_axis) / norm(cross(N, x_axis))

    # Project point onto the plane
    v = pos - v0
    x = dot(v, x_axis)
    y = dot(v, y_axis)

    return PVector2D(x, y)
end

"""
$(TYPEDSIGNATURES)
Ray casting method for Polygon3D - project to plane then use 2D ray casting.
"""
function is_inbound_ray_casting(pos::PVector, polygon::Polygon3D)
    # First check if the point is on the plane
    vertices = polygon.vertices
    N = normal(polygon)
    v0 = vertices[1]
    
    # Distance from point to plane
    dist = abs(dot(pos - v0, N))
    threshold = unit(pos.x)^2 * 1e-6
    if dist > threshold
        return false
    end
    
    # Project to 2D and use ray casting
    pos_2d = project_to_plane(pos, polygon)
    vertices_2d = [project_to_plane(v, polygon) for v in vertices]
    polygon_2d = Polygon2D(vertices_2d)
    
    return is_inbound_ray_casting(pos_2d, polygon_2d)
end

"""
$(TYPEDSIGNATURES)
For convex polygons, use cross product method to check whether the point is inside the polygon.
Points lying on an edge or vertex produce a zero cross product for the
adjacent edges and are therefore considered "in" the polygon.
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
        s = sign(dot(cross(edge, vp), N))
        if s == 0
            # ``pos`` lies on the edge ``v1 - v2``; the sign test is
            # uninformative here, so just skip it.
            continue
        elseif cross_sign == 0
            cross_sign = s
        elseif s != cross_sign
            return false
        end
    end
    return true
end

function is_inbound(pos::PVector2D, polygon::Polygon2D)
    if isconvex(polygon)
        return is_inbound_cross_product(pos, polygon)
    else
        return is_inbound_ray_casting(pos, polygon)
    end
end

function is_inbound(pos::PVector, polygon::Polygon3D)
    # First, the point must lie on (or very close to) the polygon's plane,
    # otherwise the 2D-projection based checks below would silently report
    # the point as being inside a polygon that it is actually not on.
    vertices = polygon.vertices
    v0 = vertices[1]
    N = normal(polygon)
    dist_to_plane = abs(dot(pos - v0, N))
    threshold = unit(pos.x)^2 * 1e-6
    if dist_to_plane > threshold
        return false
    end

    if isconvex(polygon)
        return is_inbound_cross_product(pos, polygon)
    else
        return is_inbound_ray_casting(pos, polygon)
    end
end

### Some common polygons

"""
$(TYPEDSIGNATURES)
Create a 2D rectangle polygon centered at origin with given width and height.
"""
function polygon_rect(width::T, height::U) where {T, U}
    w = convert(promote_type(T, U), width) / 2
    h = convert(promote_type(T, U), height) / 2
    return Polygon2D([
        PVector2D(w, h),
        PVector2D(w, -h),
        PVector2D(-w, -h),
        PVector2D(-w, h)
    ])
end

"""
$(TYPEDSIGNATURES)
Create a regular n-sided polygon with given circumradius.
"""
function polygon_regular(n::Int, radius::T) where T
    if n < 3
        error("A polygon must have at least 3 sides")
    end
    vertices = PVector2D{typeof(radius)}[]
    for i in 0:n-1
        angle = 2π * i / n
        push!(vertices, PVector2D(radius * cos(angle), radius * sin(angle)))
    end
    return Polygon2D(vertices)
end