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

#orient_exact(t::AbstractTriangle) = orient_exact(t.a, t.b, t.c)

area(t::AbstractTriangle2D) = 0.5 * abs(getz(orient(t)))

area(t::AbstractTriangle3D) = 0.5 * norm(orient(t))

centroid(t::AbstractTriangle) = (t.a + t.b + t.c) / 3.0

### Refer to Meshkit project

circumcenter(t::AbstractTriangle) = circumcenter(t.a, t.b, t.c)

#circumcenter_exact(t::AbstractTriangle) = circumcenter_exact(t.a, t.b, t.c)

incircle(t::AbstractTriangle2D, p::AbstractPoint2D) = incircle(t.a, t.b, t.c, p)

#incircle_exact(t::AbstractTriangle2D, p::AbstractPoint2D) = incircle_exact(t.a, t.b, t.c, p)

function orientation(t::AbstractTriangle2D)
    x = ustrip(getz(orient(t)))
    if x < 0
        return NegativelyOriented()
    elseif x > 0
        return PositivelyOriented()
    end
    return UnOriented()
end
#=
function orientation_exact(t::AbstractTriangle2D)
    x = ustrip(getz(orient_exact(t)))
    if x < 0
        return NegativelyOriented()
    elseif x > 0
        return PositivelyOriented()
    end
    return UnOriented()
end=#