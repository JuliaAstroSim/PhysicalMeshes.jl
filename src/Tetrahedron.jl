struct Tetrahedron{T<:Number} <: AbstractTetrahedron{T}
    a::PVector{T}
    b::PVector{T}
    c::PVector{T}
    d::PVector{T}
end

centroid(t::AbstractTetrahedron) = (t.a + t.b + t.c + t.d) / 4.0

circumcenter(t::AbstractTetrahedron) = circumcenter(t.a, t.b, t.c, t.d)

insphere(t::AbstractTetrahedron, p::AbstractPoint3D) = insphere(t.a, t.b, t.c, t.d, p)

orient(t::AbstractTetrahedron) = orient(t.a, t.b, t.c, t.d)

volume(t::AbstractTetrahedron) = norm(orient(t)) / 6.0

function orientation(t::AbstractTetrahedron)
    x = ustrip(orient(t))
    if x < 0
        return NegativelyOriented()
    elseif x > 0
        return PositivelyOriented()
    end
    return UnOriented()
end