function incircle(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D, d::AbstractPoint2D)
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

    result = ustrip(alift * bcdet + blift * cadet + clift * abdet)

    if result < 0
        return Inner()
    elseif result > 0
        return Outter()
    end

    return OnEdge()
end

function insphere(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D, e::AbstractPoint3D)
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

    result = ustrip((dlift * abc - clift * dab) + (blift * cda - alift * bcd))

    if result < 0
        return Inner()
    elseif result > 0
        return Outter()
    end

    return OnEdge()
end