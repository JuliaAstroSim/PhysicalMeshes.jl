# Refer to Meshkit project

"""
function orient(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)

    Computes the orient.
    The result is also a rough approximation of twice the signed area.
"""
function orient(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)
    acx = a.x - c.x
    bcx = b.x - c.x
    acy = a.y - c.y
    bcy = b.y - c.y
    z = acx * bcy - acy * bcx
    return PVector(zero(z), zero(z), z)
end

"""
function orient(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D)

    Computes the orient.
    The result is also a rough approximation of twice the signed area.
"""
function orient(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D)
    acx = a.x - c.x
    bcx = b.x - c.x
    acy = a.y - c.y
    bcy = b.y - c.y
    acz = a.z - c.z
    bcz = b.z - c.z
    return PVector(
        acy * bcz - acz * bcy,
        acx * bcz - acz * bcx,
        acx * bcy - acy * bcx
    )
end

"""
function orient(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)

    Computes the exact orient..
    The result is also a rough approximation of six times the signed volume.
"""
function orient(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    adx = a.x - d.x
    bdx = b.x - d.x
    cdx = c.x - d.x
    ady = a.y - d.y
    bdy = b.y - d.y
    cdy = c.y - d.y
    adz = a.z - d.z
    bdz = b.z - d.z
    cdz = c.z - d.z
    return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) + cdx * (ady * bdz - adz * bdy)
end