"""
    struct MeshConfig

- `mode`: the way of sampling physical properties
  - `CellMode`: properties are located in cell centers. `Nx × Ny × Nz` data points in total.
  - `VertexMode`: properties are located on grid points. `(Nx + 1) × (Ny + 1) × (Nz + 1)` data points in total.
- `Nx`, `Ny`, `Nz`: number of cells. In `VertexMode`, there are N+1 vertices in each direction
"""
struct MeshConfig{I,LEN,VI,V}
    mode::MeshMode
    assignment::MeshAssignment
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

    # Vector info
    Δ::V
    Min::V
    Max::V
    N::VI
end

function MeshConfig(units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
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

    Δ = SVector(Δx, Δy, Δz)
    Min = SVector(xMin, yMin, zMin)
    Max = SVector(xMax, yMax, zMax)
    N = SVector(Nx, Ny, Nz)
    return MeshConfig(mode,assignment,Nx,Ny,Nz,xMin,xMax,yMin,yMax,zMin,zMax,Δx,Δy,Δz, Δ,Min,Max,N)
end

function MeshConfig(e::Extent, units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
)
    xMin = e.xMin
    yMin = e.yMin
    zMin = e.zMin
    xMax = e.xMax
    yMax = e.yMax
    zMax = e.zMax
    Δx = (xMax-xMin)/Nx
    Δy = (yMax-yMin)/Ny
    Δz = (zMax-zMin)/Nz

    Δ = SVector(Δx, Δy, Δz)
    Min = SVector(xMin, yMin, zMin)
    Max = SVector(xMax, yMax, zMax)
    N = SVector(Nx, Ny, Nz)
    return MeshConfig(mode,assignment,Nx,Ny,Nz,xMin,xMax,yMin,yMax,zMin,zMax,Δx,Δy,Δz, Δ,Min,Max,N)
end

struct MeshCartesianStatic{I, LEN, POS, VEL, ACC, _e, RHO, PHI, _B, _E, _U, _F, _G, _H, _J}
    config::MeshConfig{I, LEN}
    pos::POS
    vel::VEL
    acc::ACC
    e::_e        # energy
    rho::RHO    # density
    phi::PHI    # potential

    B::_B       # magnetic field
    E::_E       # eletric field

    # CFD
    U::_U
    F::_F
    G::_G
    H::_H
    J::_J
end

function __MeshCartesianStatic(config::MeshConfig, ::VertexMode, units = nothing)
    x = collect(config.xMin:config.Δx:config.xMax)
    y = collect(config.yMin:config.Δy:config.yMax)
    z = collect(config.zMin:config.Δz:config.zMax)

    iter = Iterators.product(x,y,z)

    uLength = getuLength(units)

    zv = ZeroValue(units)

    pos = StructArray(PVector(p..., uLength) for p in iter)
    vel = StructArray(zv.vel for p in iter)
    acc = StructArray(zv.acc for p in iter)
    e = MArray{Tuple{(config.N.+1)...}}([zv.potpermass for p in iter])
    rho = MArray{Tuple{(config.N.+1)...}}([zv.density for p in iter])
    phi = MArray{Tuple{(config.N.+1)...}}([zv.potpermass for p in iter])

    return MeshCartesianStatic(
        config,
        pos, vel, acc, e, rho, phi,
        nothing, nothing,
        nothing, nothing, nothing, nothing, nothing,
    )
end

function __MeshCartesianStatic(config::MeshConfig, ::CellMode, units = nothing)

end

function MeshCartesianStatic(config::MeshConfig, units = nothing)
    return __MeshCartesianStatic(config, config.mode, units)
end

function MeshCartesianStatic(particles::StructArray, units = nothing;
    mode = VertexMode(),    
    assignment = CIC(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
    assign = true,
    kw...
)
    e = extent(particles)
    config = MeshConfig(e, units; Nx, Ny, Nz, mode)
    mesh = __MeshCartesianStatic(config, mode, units)

    if assign
        assignmesh(particles, mesh, assignment)
    end
    return mesh
end

function MeshCartesianStatic(particles::Array, units = nothing; kw...)
    return MeshCartesianStatic(StructArray(particles))
end