"""
    struct MeshConfig

- `mode`: the way of sampling physical properties
  - `CellMode`: properties are located in cell centers. `Nx × Ny × Nz` data points in total.
  - `VertexMode`: properties are located on grid points. `(Nx + 1) × (Ny + 1) × (Nz + 1)` data points in total.
- `Nx`, `Ny`, `Nz`: number of cells. In `VertexMode`, there are N+1 vertices in each direction
"""
struct MeshConfig{I,LEN, MODE}
    mode::MODE
    Nx::I
    Ny::I
    Nz::I

    xMin::LEN
    xMax::LEN
    yMin::LEN
    yMax::LEN
    zMin::LEN
    zMax::LEN
    Δx::LEN
    Δy::LEN
    Δz::LEN
end

function MeshConfig(units = nothing;
        mode = VertexMode(),
        Nx = 10,
        Ny = 10,
        Nz = 10,
        xMin = isnothing(units) ? -1.0 : -1.0 * units[1],
        xMax = isnothing(units) ? +1.0 : +1.0 * units[1],
        yMin = isnothing(units) ? -1.0 : -1.0 * units[1],
        yMax = isnothing(units) ? +1.0 : +1.0 * units[1],
        zMin = isnothing(units) ? -1.0 : -1.0 * units[1],
        zMax = isnothing(units) ? +1.0 : +1.0 * units[1],
    )
        Δx = (xMax-xMin)/Nx
        Δy = (yMax-yMin)/Ny
        Δz = (zMax-zMin)/Nz
    return MeshConfig(mode,Nx,Ny,Nz,xMin,xMax,yMin,yMax,zMin,zMax,Δx,Δy,Δz)
end

struct MeshCartesianStatic{I, LEN, POS, VEL, ACC, E, RHO, PHI, _U, _F, _G, _H, _J}
    config::MeshConfig{I, LEN}
    pos::POS
    vel::VEL
    acc::ACC
    e::E        # energy
    rho::RHO    # density
    phi::PHI    # potential

    # CFD
    U::_U
    F::_F
    G::_G
    H::_H
    J::_J
end

function MeshCartesianStatic(mode::VertexMode, units = nothing;
        config = MeshConfig(units; mode),
    )
    x = collect(config.xMin:config.Δx:config.xMax)
    y = collect(config.yMin:config.Δy:config.yMax)
    z = collect(config.zMin:config.Δz:config.zMax)

    iter = Iterators.product(x,y,z)

    uLength = getuLength(units)

    zv = ZeroValue(units)

    pos = StructArray(PVector(p..., uLength) for p in iter)
    vel = StructArray(zv.vel for p in iter)
    acc = StructArray(zv.acc for p in iter)
    e = StructArray(zv.potpermass for p in iter)
    rho = StructArray(zv.density for p in iter)
    phi = StructArray(zv.potpermass for p in iter)

    return MeshCartesianStatic(
        config,
        pos, vel, acc, e, rho, phi,
        nothing, nothing, nothing, nothing, nothing,
    )
end