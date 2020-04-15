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

#function insphere(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
#    
#end