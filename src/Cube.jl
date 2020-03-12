struct Cube2D{T<:Number} <: AbstractCube2D{T}
    top::PVector2D{T}
    below::PVector2D{T}
end

struct Cube{T<:Number} <: AbstractCube3D{T}
    top::PVector{T}
    below::PVector{T}
end

Cube(top::PVector2D{T}, below::PVector2D{T}) where T<:Number = Cube2D(top, below)

#! Area and Volume is signed

function area(c::AbstractCube2D)
    p = c.top - c.below
    return p.x * p.y
end

function volume(c::AbstractCube3D)
    p = c.top - c.below
    return p.x * p.y * p.z
end