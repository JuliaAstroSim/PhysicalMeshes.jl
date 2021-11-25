"""
    struct Tetrahedron{T<:Number} <: AbstractTetrahedron{T}

## Fields
- a::PVector{T}
- b::PVector{T}
- c::PVector{T}
- d::PVector{T}

## Examples
```julia
julia> t = Tetrahedron(
           PVector(1.0, 0.0, 1.0),
           PVector(1.0, 1.0, 0.0),
           PVector(0.0, 1.0, 1.0),
           PVector(1.0, 1.0, 1.0)
       )
Tetrahedron{Float64}(PVector{Float64}(1.0, 0.0, 1.0), PVector{Float64}(1.0, 1.0, 0.0), PVector{Float64}(0.0, 1.0, 1.0), PVector{Float64}(1.0, 1.0, 1.0))

julia> centroid(t)
PVector{Float64}(0.75, 0.75, 0.75)

julia> circumcenter(t)
PVector{Float64}(0.5, 0.5, 0.5)

julia> insphere(t, PVector(0.5, 0.5, 0.5))
Interior()

julia> insphere(t, PVector())
OnEdge()

julia> insphere(t, PVector(2.0, 0.0, 0.0))
Exterior()

julia> PhysicalMeshes.volume(t)
0.16666666666666666

julia> orientation(t)
PositivelyOriented()
```
"""
struct Tetrahedron{T<:Number} <: AbstractTetrahedron{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
    d::PVector{T}
end

"""
    centroid(t::AbstractTetrahedron)

Averaged position center of the tetrahedron
"""
centroid(t::AbstractTetrahedron) = (t.a + t.b + t.c + t.d) / 4.0

"""
    circumcenter(t::AbstractTetrahedron)

Return the circumsphere center of tetrahedron
"""
circumcenter(t::AbstractTetrahedron) = circumcenter(t.a, t.b, t.c, t.d)

#circumcenter_exact(t::AbstractTetrahedron) = circumcenter_exact(t.a, t.b, t.c, t.d)

"""
insphere(t::AbstractTetrahedron, p::AbstractPoint3D)

Test whether a point locates inside the circumsphere center of tetrahedron. Return a trait (`Interior`, `Exterior`, or `OnEdge`)
"""
insphere(t::AbstractTetrahedron, p::AbstractPoint3D) = insphere(t.a, t.b, t.c, t.d, p)

#insphere_exact(t::AbstractTetrahedron, p::AbstractPoint3D) = insphere_exact(t.a, t.b, t.c, t.d, p)

"""
    orient(t::AbstractTetrahedron)

Return a positive number if positively oriented.
Return a negative number if negatively oriented.
Return zero if not oriented (indicating the tetrahedron has zero volume)
"""
orient(t::AbstractTetrahedron) = orient(t.a, t.b, t.c, t.d)

#orient_exact(t::AbstractTetrahedron) = orient_exact(t.a, t.b, t.c, t.d)

"""
    volume(t::AbstractTetrahedron)

Unsigned volume of tetrahedron
"""
volume(t::AbstractTetrahedron) = norm(orient(t)) / 6.0

"""
    orientation(t::AbstractTetrahedron)

Orientation of tetrahedron. Return a trait (`NegativelyOriented`, `PositivelyOriented`, or `UnOriented`)
"""
function orientation(t::AbstractTetrahedron)
    x = orient(t)
    if x < 0
        return NegativelyOriented()
    elseif x > 0
        return PositivelyOriented()
    end
    return UnOriented()
end
#=
function orientation_exact(t::AbstractTetrahedron)
    x = orient_exact(t)
    if x < 0
        return NegativelyOriented()
    elseif x > 0
        return PositivelyOriented()
    end
    return UnOriented()
end=#