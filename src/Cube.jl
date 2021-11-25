"""
    struct Cube2D{T<:Number} <: AbstractCube2D{T}

## Fields
- top::PVector2D{T}
- below::PVector2D{T}

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
    struct Cube{T<:Number} <: AbstractCube3D{T}

## Fields
- top::PVector{T}
- below::PVector{T}

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

#! Area and Volume is signed

"""
    area(c::AbstractCube2D)

Signed area of 2D cube
"""
function area(c::AbstractCube2D)
    p = c.top - c.below
    return p.x * p.y
end

"""
    volume(c::AbstractCube3D)

Signed volume of 3D cube
"""
function volume(c::AbstractCube3D)
    p = c.top - c.below
    return p.x * p.y * p.z
end