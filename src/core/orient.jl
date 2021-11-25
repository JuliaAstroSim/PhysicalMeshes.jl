# Refer to Meshkit project

function orient_kernel(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)
    acx = a.x - c.x
    bcx = b.x - c.x
    acy = a.y - c.y
    bcy = b.y - c.y
    return acx * bcy - acy * bcx
end

"""
function orient(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)

    Computes the orient.
    The result is also a rough approximation of twice the signed area.
"""
function orient(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)
    z = orient_kernel(a, b, c)
    return PVector(zero(z), zero(z), z)
end
#=
function orient_exact(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)
    z = floatnumber(orient_kernel(decimal(a), decimal(b), decimal(c)))
    return PVector(zero(z), zero(z), z)
end=#

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
#=
function orient_exact(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D)
    ea = decimal(a)
    eb = decimal(b)
    ec = decimal(c)
    acx = ea.x - ec.x
    bcx = eb.x - ec.x
    acy = ea.y - ec.y
    bcy = eb.y - ec.y
    acz = ea.z - ec.z
    bcz = eb.z - ec.z
    return floatnumber(PVector(
        acy * bcz - acz * bcy,
        acx * bcz - acz * bcx,
        acx * bcy - acy * bcx
    ))
end
=#
function orient_kernel(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    ax = a.x - d.x
    bx = b.x - d.x
    cx = c.x - d.x
    ay = a.y - d.y
    by = b.y - d.y
    cy = c.y - d.y
    az = a.z - d.z
    bz = b.z - d.z
    cz = c.z - d.z
    return ustrip(ax * (bz * cy - by * cz) + bx * (cz * ay - cy * az) + cx * (az * by - ay * bz))
end

"""
function orient(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)

    Computes the orient.
    The result is also a rough approximation of six times the signed volume.
"""
function orient(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    return orient_kernel(a, b, c, d)
end
#=
function orient_exact(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    return orient_kernel(decimal(a), decimal(b), decimal(c), decimal(d))
end=#