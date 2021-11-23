function incircle_kernel(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D, d::AbstractPoint2D)
    adx = a.x - d.x
    ady = a.y - d.y
    bdx = b.x - d.x
    bdy = b.y - d.y
    cdx = c.x - d.x
    cdy = c.y - d.y

    abdet = adx * bdy - bdx * ady
    bcdet = bdx * cdy - cdx * bdy
    cadet = cdx * ady - adx * cdy
    alift = adx * adx + ady * ady
    blift = bdx * bdx + bdy * bdy
    clift = cdx * cdx + cdy * cdy

    return alift * bcdet + blift * cadet + clift * abdet
end

function incircle(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D, d::AbstractPoint2D)
    result = ustrip(incircle_kernel(a, b, c, d))

    if result < 0
        return Interior()
    elseif result > 0
        return Exterior()
    end

    return OnEdge()
end
#=
function insphere_exact(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D, d::AbstractPoint2D)
    result = ustrip(floatnumber(incircle_kernel(decimal(a), decimal(b), decimal(c), decimal(d))))

    if result < 0
        return Interior()
    elseif result > 0
        return Exterior()
    end

    return OnEdge()
end=#

function insphere_kernel(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D, e::AbstractPoint3D)
    aex = a.x - e.x
    bex = b.x - e.x
    cex = c.x - e.x
    dex = d.x - e.x
    aey = a.y - e.y
    bey = b.y - e.y
    cey = c.y - e.y
    dey = d.y - e.y
    aez = a.z - e.z
    bez = b.z - e.z
    cez = c.z - e.z
    dez = d.z - e.z

    ab = aex * bey - bex * aey
    bc = bex * cey - cex * bey
    cd = cex * dey - dex * cey
    da = dex * aey - aex * dey

    ac = aex * cey - cex * aey
    bd = bex * dey - dex * bey

    abc = aez * bc - bez * ac + cez * ab
    bcd = bez * cd - cez * bd + dez * bc
    cda = cez * da + dez * ac + aez * cd
    dab = dez * ab + aez * bd + bez * da

    alift = aex * aex + aey * aey + aez * aez
    blift = bex * bex + bey * bey + bez * bez
    clift = cex * cex + cey * cey + cez * cez
    dlift = dex * dex + dey * dey + dez * dez

    return (dlift * abc - clift * dab) + (blift * cda - alift * bcd)
end

function insphere(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D, e::AbstractPoint3D)
    result = ustrip(insphere_kernel(a, b, c, d, e))

    if result < 0
        return Interior()
    elseif result > 0
        return Exterior()
    end

    return OnEdge()
end
#=
function insphere_exact(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D, e::AbstractPoint3D)
    result = ustrip(floatnumber(insphere_kernel(decimal(a), decimal(b), decimal(c), decimal(d), decimal(e))))

    if result < 0
        return Interior()
    elseif result > 0
        return Exterior()
    end

    return OnEdge()
end=#