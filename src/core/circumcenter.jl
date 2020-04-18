function circumcenter(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)
    xba = b.x - a.x
    yba = b.y - a.y
    xca = c.x - a.x
    yca = c.y - a.y
    balength = xba * xba + yba * yba
    calength = xca * xca + yca * yca

    denominator = 0.5 / (xba * yca - yba * xca)

    return PVector2D(
        a.x + (yca * balength - yba * calength) * denominator,
        a.y + (xba * calength - xca * balength) * denominator
    )
end

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

    denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc)

    return PVector(
        a.x + ((balength * yca - calength * yba) * zcrossbc - (balength * zca - calength * zba) * ycrossbc) * denominator,
        a.y + ((balength * zca - calength * zba) * xcrossbc - (balength * xca - calength * xba) * zcrossbc) * denominator,
        a.z + ((balength * xca - calength * xba) * ycrossbc - (balength * yca - calength * yba) * xcrossbc) * denominator
    )   
end

function circumcenter_exact(a::AbstractPoint, b::AbstractPoint, c::AbstractPoint)
    return number(circumcenter(decimal(a), decimal(b), decimal(c)))
end

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
    denominator = 0.5 / (xba * xcrosscd + yba * ycrosscd + zba * zcrosscd)
    return PVector(
        a.x + (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) * denominator,
        a.y + (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) * denominator,
        a.z + (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) * denominator
    )
end

function circumcenter_exact(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    return number(circumcenter(decimal(a), decimal(b), decimal(c), decimal(d)))
end