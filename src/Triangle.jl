struct Triangle2D{T<:Number} <: AbstractTriangle2D{T}
    a::PVector2D{T}
    b::PVector2D{T}
    c::PVector2D{T}
end

struct Triangle{T<:Number} <: AbstractTriangle3D{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
end

# https://en.wikipedia.org/wiki/Triangle

len(t::AbstractTriangle) = norm(t.a - t.b) + norm(t.a - t.c) + norm(t.b - t.c)

orient(t::AbstractTriangle) = orient(t.a, t.b, t.c)

area(t::AbstractTriangle2D) = 0.5 * abs(getz(orient(t)))

area(t::AbstractTriangle3D) = 0.5 * norm(orient(t))

orientation(t::AbstractTriangle2D) = ustrip(getz(orient(t))) >= 0 ? PositivelyOriented() : NegativelyOriented()



centroid(t::AbstractTriangle) = (t.a + t.b + t.c) / 3.0

### Refer to Meshkit project

circumcenter(t::AbstractTriangle) = circumcenter(t.a, t.b, t.c)

#circumcenter_exact(t::AbstractTriangle2D)

#circumcenter_exact(t::AbstractTriangle3D)


incircle(t::AbstractTriangle2D, p::AbstractPoint2D) = incircle(t.a, t.b, t.c, p)

#incircle_exact()

#incircle_exact()
