"""
    struct Triangle2D{T<:Number} <: AbstractTriangle2D{T}

## Fields
- a::PVector2D{T}
- b::PVector2D{T}
- c::PVector2D{T}

## Examples
```julia
julia> t = Triangle2D(PVector2D(4.0, 3.0), PVector2D(4.0, 0.0), PVector2D(0.0, 3.0))
Triangle2D{Float64}(PVector2D{Float64}(4.0, 3.0), PVector2D{Float64}(4.0, 0.0), PVector2D{Float64}(0.0, 3.0))

julia> orientation(t)
NegativelyOriented()

julia> len(t)
12.0

julia> PhysicalMeshes.area(t)
6.0

julia> centroid(t)
PVector2D{Float64}(2.6666666666666665, 2.0)

julia> circumcenter(t)
PVector2D{Float64}(2.0, 1.5)

julia> incircle(t, PVector2D(3.0, 2.0))
Interior()

julia> incircle(t, PVector2D())
OnEdge()

julia> incircle(t, PVector2D(6.0, 0.0))
Exterior()
```
"""
struct Triangle2D{T<:Number} <: AbstractTriangle2D{T}
    a::PVector2D{T}
    b::PVector2D{T}
    c::PVector2D{T}
end

"""
    struct Triangle{T<:Number} <: AbstractTriangle3D{T}

## Fields
- a::PVector{T}
- b::PVector{T}
- c::PVector{T}

## Examples
```julia
julia> t = Triangle(PVector(0.0, 0.0, 0.0, u"m"), PVector(3.0, 4.0, 0.0, u"m"), PVector(3.0, 4.0, 12.0, u"m"))
Triangle{Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}(PVector(0.0 m, 0.0 m, 0.0 m), PVector(3.0 m, 4.0 m, 0.0 m), PVector(3.0 m, 4.0 m, 12.0 m))

julia> len(t)
30.0 m

julia> PhysicalMeshes.area(t)
30.0 m^2

julia> centroid(t)
PVector(2.0 m, 2.6666666666666665 m, 4.0 m)

julia> circumcenter(t)
PVector(1.5 m, 2.0 m, 6.0 m)
```
"""
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