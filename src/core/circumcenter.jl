"""
$(TYPEDSIGNATURES)
Compute the circumcenter of three 2D points.

Throws an `ArgumentError` if the three points are collinear (no unique
circumcenter exists).
"""
function circumcenter(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)
    xba = b.x - a.x
    yba = b.y - a.y
    xca = c.x - a.x
    yca = c.y - a.y
    balength = xba * xba + yba * yba
    calength = xca * xca + yca * yca

    denom = xba * yca - yba * xca
    if isapprox(denom, 0; atol = eps(typeof(denom)) * 100)
        throw(ArgumentError("Collinear points: no unique circumcenter"))
    end
    denominator = 0.5 / denom

    return PVector2D(
        a.x + (yca * balength - yba * calength) * denominator,
        a.y + (xba * calength - xca * balength) * denominator
    )
end

"""
$(TYPEDSIGNATURES)
Compute the circumcenter of three 3D points.

Throws an `ArgumentError` if the three points are collinear (no unique
circumcenter exists).
"""
function circumcenter(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D)
    xba = b.x - a.x
    yba = b.y - a.y
    zba = b.z - a.z
    xca = c.x - a.x
    yca = c.y - a.y
    zca = c.z - a.z
    balength = xba * xba + yba * yba + zba * zba
    calength = xca * xca + yca * yca + zca * zca
    xcrossbc = yba * zca - yca * zba
    ycrossbc = zba * xca - zca * xba
    zcrossbc = xba * yca - xca * yba

    denom = xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc
    if isapprox(denom, 0; atol = eps(typeof(denom)) * 100)
        throw(ArgumentError("Collinear points: no unique circumcenter"))
    end
    denominator = 0.5 / denom

    return PVector(
        a.x + ((balength * yca - calength * yba) * zcrossbc - (balength * zca - calength * zba) * ycrossbc) * denominator,
        a.y + ((balength * zca - calength * zba) * xcrossbc - (balength * xca - calength * xba) * zcrossbc) * denominator,
        a.z + ((balength * xca - calength * xba) * ycrossbc - (balength * yca - calength * yba) * xcrossbc) * denominator
    )
end

function circumcenter_exact(a::AbstractPoint, b::AbstractPoint, c::AbstractPoint)
    return floatnumber(circumcenter(decimal(a), decimal(b), decimal(c)))
end

"""
$(TYPEDSIGNATURES)
Compute the circumcenter of four 3D points (the center of the circumscribed
sphere of a tetrahedron).

Throws an `ArgumentError` if the four points are coplanar (no unique
circumcenter exists).
"""
function circumcenter(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    xba = b.x - a.x
    yba = b.y - a.y
    zba = b.z - a.z
    xca = c.x - a.x
    yca = c.y - a.y
    zca = c.z - a.z
    xda = d.x - a.x
    yda = d.y - a.y
    zda = d.z - a.z
    balength = xba * xba + yba * yba + zba * zba
    calength = xca * xca + yca * yca + zca * zca
    dalength = xda * xda + yda * yda + zda * zda
    xcrosscd = yca * zda - yda * zca
    ycrosscd = zca * xda - zda * xca
    zcrosscd = xca * yda - xda * yca
    xcrossdb = yda * zba - yba * zda
    ycrossdb = zda * xba - zba * xda
    zcrossdb = xda * yba - xba * yda
    xcrossbc = yba * zca - yca * zba
    ycrossbc = zba * xca - zca * xba
    zcrossbc = xba * yca - xca * yba
    denom = xba * xcrosscd + yba * ycrosscd + zba * zcrosscd
    if isapprox(denom, 0; atol = eps(typeof(denom)) * 100)
        throw(ArgumentError("Coplanar points: no unique circumcenter"))
    end
    denominator = 0.5 / denom
    return PVector(
        a.x + (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) * denominator,
        a.y + (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) * denominator,
        a.z + (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) * denominator
    )
end

function circumcenter_exact(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    return floatnumber(circumcenter(decimal(a), decimal(b), decimal(c), decimal(d)))
end
