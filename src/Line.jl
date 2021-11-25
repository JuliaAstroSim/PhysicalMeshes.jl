"""
    struct Line2D{T<:Number} <: AbstractLine2D{T}

## Fields
- a::PVector2D{T}
- b::PVector2D{T}

## Examples
```julia
julia> line = Line(PVector2D(), PVector2D(1.0, 0.0))
Line2D{Float64}(PVector2D{Float64}(0.0, 0.0), PVector2D{Float64}(1.0, 0.0))

julia> len(line)
1.0

julia> midpoint(line)
PVector2D{Float64}(0.5, 0.0)

julia> line + PVector(1.0,1.0)
Line2D{Float64}(PVector2D{Float64}(1.0, 1.0), PVector2D{Float64}(2.0, 1.0))

julia> line * 2
Line2D{Float64}(PVector2D{Float64}(0.0, 0.0), PVector2D{Float64}(2.0, 0.0))
```
"""
struct Line2D{T<:Number} <: AbstractLine2D{T}
    a::PVector2D{T}
    b::PVector2D{T}
end

@inline +(m::Line2D, n::Line2D) = Line2D(m.a + n.a, m.b + n.b)

@inline -(m::Line2D, n::Line2D) = Line2D(m.a - n.a, m.b - n.b)

@inline *(m::Line2D, h::Number) = Line2D(m.a * h, m.b * h)
@inline *(h::Number, m::Line2D) = Line2D(m.a * h, m.b * h)

#! It would be confusing to add a point with a line
@inline +(m::Line2D, p::PVector2D) = Line2D(m.a + p, m.b + p)
@inline -(m::Line2D, p::PVector2D) = Line2D(m.a - p, m.b - p)
@inline /(m::Line2D, h::Number) = Line2D(m.a / h, m.b / h)

"""
    struct Line{T<:Number} <: AbstractLine3D{T}

## Fields
- a::PVector{T}
- b::PVector{T}

## Examples
```julia
julia> line = Line(PVector(u"m"), PVector(1.0, 0.0, 0.0, u"m"))
Line{Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}(PVector(0.0 m, 0.0 m, 0.0 m), PVector(1.0 m, 0.0 m, 0.0 m))

julia> len(line)
1.0 m

julia> midpoint(line)
PVector(0.5 m, 0.0 m, 0.0 m)

julia> line + PVector(1.0,2.0,3.0,u"m")
Line{Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}(PVector(1.0 m, 2.0 m, 3.0 m), PVector(2.0 m, 2.0 m, 3.0 m))

julia> line * 3
Line{Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}(PVector(0.0 m, 0.0 m, 0.0 m), PVector(3.0 m, 0.0 m, 0.0 m))
```
"""
struct Line{T<:Number} <: AbstractLine3D{T}
    a::PVector{T}
    b::PVector{T}
end

Line(a::AbstractPoint2D{T}, b::AbstractPoint2D{T}) where T<:Number = Line2D(a,b)

@inline +(m::Line, n::Line) = Line(m.a + n.a, m.b + n.b)

@inline -(m::Line, n::Line) = Line(m.a - n.a, m.b - n.b)

@inline *(m::Line, h::Number) = Line(m.a * h, m.b * h)
@inline *(h::Number, m::Line) = Line(m.a * h, m.b * h)

#! It would be confusing to add a point with a line
@inline +(m::Line, p::PVector) = Line(m.a + p, m.b + p)
@inline -(m::Line, p::PVector) = Line(m.a - p, m.b - p)
@inline /(m::Line, h::Number) = Line(m.a / h, m.b / h)

# General functions

"""
    len(line::AbstractLine)

Length of line
"""
len(line::AbstractLine) = norm(line.a - line.b)

"""
    centroid(line::AbstractLine)

Average position of line vertices.
"""
centroid(line::AbstractLine) = (line.a + line.b) * 0.5

"""
    center(line::AbstractLine) = centroid(line)

Average position of line vertices.
"""
center(line::AbstractLine) = centroid(line)

"""
    midpoint(line::AbstractLine) = centroid(line)

Average position of line vertices.
"""
midpoint(line::AbstractLine) = centroid(line)