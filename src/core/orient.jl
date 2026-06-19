# Refer to Meshkit project

"""
$(TYPEDSIGNATURES)
Kernel of the 2D orientation test. Returns twice the signed area of the
triangle (a, b, c), or zero if the three points are (numerically) collinear.
"""
@inline function orient_kernel(a::AbstractPoint2D, b::AbstractPoint2D, c::AbstractPoint2D)
    acx = a.x - c.x
    bcx = b.x - c.x
    acy = a.y - c.y
    bcy = b.y - c.y
    return acx * bcy - acy * bcx
end

"""
$(TYPEDSIGNATURES)
Computes the 2D orientation. The result is also a rough approximation of
twice the signed area. Returns zero when the three points are collinear.
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
$(TYPEDSIGNATURES)
Computes the 3D orientation of the triangle (a, b, c) and returns the
normal vector whose magnitude is twice the signed area. For collinear
points the returned vector has zero magnitude.
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

"""
$(TYPEDSIGNATURES)
Kernel of the 3D orient test. Returns six times the signed volume of the
tetrahedron (a, b, c, d). The result is exactly zero when the four points
are coplanar.
"""
@inline function orient_kernel(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    ax = a.x - d.x
    bx = b.x - d.x
    cx = c.x - d.x
    ay = a.y - d.y
    by = b.y - d.y
    cy = c.y - d.y
    az = a.z - d.z
    bz = b.z - d.z
    cz = c.z - d.z
    
    # Preserve units so callers like ``volume`` and ``orientation`` work
    # for both unitful and unitless inputs. The orientation trait
    # compares against zero, which works for ``Quantity`` as well.
    return ax * (bz * cy - by * cz) + bx * (cz * ay - cy * az) + cx * (az * by - ay * bz)
    # return ustrip(ax * (bz * cy - by * cz) + bx * (cz * ay - cy * az) + cx * (az * by - ay * bz))
end

"""
$(TYPEDSIGNATURES)
Computes the 3D orient. The result is a rough approximation of six times
the signed volume of the tetrahedron (a, b, c, d). The value is zero
when the four points are coplanar.
"""
@inline function orient(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    return orient_kernel(a, b, c, d)
end
#=
function orient_exact(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
    return orient_kernel(decimal(a), decimal(b), decimal(c), decimal(d))
end=#
