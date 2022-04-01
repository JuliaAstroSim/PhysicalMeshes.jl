"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct Sphere{P, L}
    center::P
    radius::L
end

"""
$(TYPEDSIGNATURES)
Return `true` if the distance of `p` from `s.center` is **smaller** than `s.radius`.
Otherwise (>=), return `false`.
"""
function interior(s::Sphere, p)
    dr = norm(p - s.center)
    if dr < s.radius
        return true
    else
        return false
    end
end

"""
$(TYPEDSIGNATURES)
Return `true` if the distance of `p` from `s.center` is **larger** than `s.radius`.
Otherwise (<=), return `false`.
"""
function exterior(s::Sphere, p)
    dr = norm(p - s.center)
    if dr > s.radius
        return true
    else
        return false
    end
end