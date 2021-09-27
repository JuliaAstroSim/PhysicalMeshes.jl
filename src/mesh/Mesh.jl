"""
    struct MeshConfig

- `mode`: the way of sampling physical properties
  - `CellMode`: properties are located in cell centers. `Nx × Ny × Nz` data points in total.
  - `VertexMode`: properties are located on grid points. `(Nx + 1) × (Ny + 1) × (Nz + 1)` data points in total.
- `Nx`, `Ny`, `Nz`: number of cells. In `VertexMode`, there are N+1 vertices in each direction
- `NB`: Number of boundary points
"""
struct MeshConfig{I,LEN,VI,V,U}
    mode::MeshMode
    assignment::MeshAssignment
    units::U
    Nx::I
    Ny::I
    Nz::I
    NB::I # Number of boundary points

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

function Base.show(io::IO, config::MeshConfig)
    print(io, 
    """
    Mesh config:
                         mode: $(config.mode)
            assignment method: $(config.assignment)
                        units: $(config.units)
              Number of Cells: $(config.N)
    Number of boundary points: $(config.NB)
                          Min: $(config.Min)
                          Max: $(config.Max)
                            Δ: $(config.Δ)
    """)
end

function MeshConfig(units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
    NB = 1,
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
    return MeshConfig(
        mode,assignment,units,
        Nx,Ny,Nz,NB,xMin,xMax,yMin,yMax,zMin,zMax,Δx,Δy,Δz,
        Δ,Min,Max,N,
    )
end

function MeshConfig(e::Extent, units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
    NB = 1,
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
    return MeshConfig(
        mode,assignment,units,
        Nx,Ny,Nz,NB,xMin,xMax,yMin,yMax,zMin,zMax,Δx,Δy,Δz,
        Δ,Min,Max,N,
    )
end

struct MeshCartesianStatic{I, LEN, VI, V, U, POS, VEL, ACC, _e, RHO, PHI, _B, _E, _U, _F, _G, _H, _J} <: AbstractMesh{U}
    config::MeshConfig{I,LEN,VI,V,U}
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
    x = collect(config.xMin - config.Δx * config.NB:config.Δx:config.xMax + config.Δx * config.NB)
    y = collect(config.yMin - config.Δy * config.NB:config.Δy:config.yMax + config.Δy * config.NB)
    z = collect(config.zMin - config.Δz * config.NB:config.Δz:config.zMax + config.Δz * config.NB)

    iter = Iterators.product(x,y,z)

    uLength = getuLength(units)

    zv = ZeroValue(units)

    pos = StructArray(PVector(p..., uLength) for p in iter)
    vel = StructArray(zv.vel for p in iter)
    acc = StructArray(zv.acc for p in iter)
    e = MArray{Tuple{(config.N.+(1+2*config.NB))...}}([zv.potpermass for p in iter])
    rho = MArray{Tuple{(config.N.+(1+2*config.NB))...}}([zv.density for p in iter])
    phi = MArray{Tuple{(config.N.+(1+2*config.NB))...}}([zv.potpermass for p in iter])

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
    NB = 1,
    assign = true,
    kw...
)
    e = extent(particles)
    config = MeshConfig(e, units; Nx, Ny, Nz, NB, mode)
    mesh = __MeshCartesianStatic(config, mode, units)

    if assign
        assignmesh(particles, mesh, assignment)
    end
    return mesh
end

function MeshCartesianStatic(particles::Array, units = nothing; kw...)
    return MeshCartesianStatic(StructArray(particles))
end