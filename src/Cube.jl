"""
$(TYPEDEF)
$(TYPEDFIELDS)

## Examples
```julia
julia> c = Cube(PVector2D(1.0, 1.0), PVector2D())
Cube2D{Float64}(PVector2D{Float64}(1.0, 1.0), PVector2D{Float64}(0.0, 0.0))

julia> PhysicalMeshes.area(c)
1.0
```
"""
struct Cube2D{T<:Number} <: AbstractCube2D{T}
    top::PVector2D{T}
    below::PVector2D{T}
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

## Examples
```julia
julia> c = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
Cube{Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}(PVector(1.0 m, 1.0 m, 1.0 m), PVector(0.0 m, 0.0 m, 0.0 m))

julia> PhysicalMeshes.volume(c)
1.0 m^3
```
"""
struct Cube{T<:Number} <: AbstractCube3D{T}
    top::PVector{T}
    below::PVector{T}
end

Cube(top::PVector2D{T}, below::PVector2D{T}) where T<:Number = Cube2D(top, below)

"""
$(TYPEDSIGNATURES)
Return `true` if `p` is inside the cube.
Otherwise (including on the edge), return `false`.
"""
function interior(c::Cube2D, p::AbstractPoint2D)
    if (c.below.x < p.x < c.top.x) && (c.below.y < p.y < c.top.y)
        return true
    else
        return false
    end
end

"""
$(TYPEDSIGNATURES)
Return `true` if `p` is inside the cube.
Otherwise (including on the edge), return `false`.
"""
function interior(c::Cube, p::AbstractPoint3D)
    if (c.below.x < p.x < c.top.x) && (c.below.y < p.y < c.top.y) && (c.below.z < p.z < c.top.z)
        return true
    else
        return false
    end
end

"""
$(TYPEDSIGNATURES)
Return `true` if `p` is outside the cube.
Otherwise (including on the edge), return `false`.
"""
function exterior(c::Cube2D, p::AbstractPoint2D)
    if (c.below.x <= p.x <= c.top.x) && (c.below.y <= p.y <= c.top.y)
        return false
    else
        return true
    end
end

"""
$(TYPEDSIGNATURES)
Return `true` if `p` is outside the cube.
Otherwise (including on the edge), return `false`.
"""
function exterior(c::Cube, p::AbstractPoint3D)
    if (c.below.x <= p.x <= c.top.x) && (c.below.y <= p.y <= c.top.y) && (c.below.z <= p.z <= c.top.z)
        return false
    else
        return true
    end
end

#! Area and Volume is signed

"""
$(TYPEDSIGNATURES)
Signed area of 2D cube
"""
function area(c::AbstractCube2D)
    p = c.top - c.below
    return p.x * p.y
end

"""
$(TYPEDSIGNATURES)
Signed volume of 3D cube
"""
function volume(c::AbstractCube3D)
    p = c.top - c.below
    return p.x * p.y * p.z
end