"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Polygon2D{T} <: AbstractPolygon2D{T}
    points
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Polygon3D{T} <: AbstractPolygon3D{T}
    points
end