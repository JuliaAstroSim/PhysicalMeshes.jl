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

len(line::AbstractLine) = norm(line.a - line.b)

centroid(line::AbstractLine) = (line.a + line.b) * 0.5

center = centroid

midpoint(line::AbstractLine) = centroid(line)