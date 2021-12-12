function mass_center(m::MeshCartesianStatic)
    rho_total = sum(m.rho)
    if iszero(rho_total)
        return mean(m.pos)
    else
        return mean(m.rho .* m.pos) / rho_total
    end
end

function total_mass(m::MeshCartesianStatic)
    return sum(m.rho .* prod(m.config.Δ))
end