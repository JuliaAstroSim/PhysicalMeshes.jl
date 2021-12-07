"""
    struct MeshConfig

## Fields
- `mode`: the way of sampling physical properties
  - `CellMode`: properties are located in cell centers. `Nx × Ny × Nz` data points in total.
  - `VertexMode`: properties are located on grid points. `(Nx + 1) × (Ny + 1) × (Nz + 1)` data points in total.
- `assignment`: mesh assignment algorithm
  - `NGP`: nearest grid point
  - `CIC`: cloud in cell
  - `TSC`: triangular shaped cloud
- `boundary`: boundary conditions
  - `Periodic`
  - `Vacuum`
  - `Dirichlet`
- `units`: support `nothing`, `uAstro`, `uSI`, `uGadget2`, `uCGS`
- `dim`: dimension of the mesh
- `NG`: Number of boundary ghost cells
- `Δ`: mesh resolution in each direction
- `Min`: minimum coordinate of each direction, not including ghost cells. Corresponds to cell coordinates in cell mode, and vertex coordinates in vertex mode.
- `Max`: maximum coordinate of each direction, not including ghost cells. Corresponds to cell coordinates in cell mode, and vertex coordinates in vertex mode.
- `N`: number of cells in each direction, not including ghost cells. In `VertexMode`, there are N+1 vertices in each direction
- `Len`: total number of cells in each direction, including ghost cells.

## Constructors
- `MeshConfig(units = nothing; kw...)`: all fields can be modified by keywords. Try construct one to see default values
- `MeshConfig(e::Extent, units = nothing; kw...)`: `Min` and `Max` are extracted from `e`. Other field can be modified by keywords

For more instructions, see the documentation: https://juliaastrosim.github.io/PhysicalMeshes.jl/dev
"""
struct MeshConfig{I,VI,V,U}
    mode::MeshMode
    assignment::MeshAssignment
    boundary::BoundaryCondition

    units::U
    dim::I
    NG::I

    # Vector info
    Δ::V
    Min::V
    Max::V
    N::VI
    Len::VI
end

function Base.show(io::IO, config::MeshConfig)
    print(io, 
    """
    Mesh config:
                          dim: $(config.dim)
                         mode: $(config.mode)
            assignment method: $(config.assignment)
           Boundary Condition: $(config.boundary)
                        units: $(config.units)
              Number of Cells: $(config.N)
       Number of ghost points: $(config.NG)
                          Min: $(config.Min)
                          Max: $(config.Max)
                            Δ: $(config.Δ)
    """)
end

function MeshConfig(units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
    boundary = Periodic(),
    Nx = 5,
    Ny = 5,
    Nz = 5,
    NG = 0,
    xMin = isnothing(units) ? -1.0 : -1.0 * units[1],
    xMax = isnothing(units) ? +1.0 : +1.0 * units[1],
    yMin = isnothing(units) ? -1.0 : -1.0 * units[1],
    yMax = isnothing(units) ? +1.0 : +1.0 * units[1],
    zMin = isnothing(units) ? -1.0 : -1.0 * units[1],
    zMax = isnothing(units) ? +1.0 : +1.0 * units[1],
    dim = 3,
)
    Δx = (xMax-xMin)/Nx
    Δy = (yMax-yMin)/Ny
    Δz = (zMax-zMin)/Nz

    Δ = SVector(Δx, Δy, Δz)
    Min = SVector(xMin, yMin, zMin)
    Max = SVector(xMax, yMax, zMax)
    N = SVector(Nx, Ny, Nz)
    Len = N .+ (2 * NG)
    return MeshConfig(
        mode,assignment,boundary,units,dim,NG,
        Δ[1:dim],Min[1:dim],Max[1:dim],N[1:dim],Len[1:dim],
    )
end

function MeshConfig(e::Extent, units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
    boundary = Periodic(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
    NG = 1,
    dim = 3,
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
    Len = N .+ (2 * NG)
    return MeshConfig(
        mode,assignment,boundary,units,dim,NG,
        Δ[1:dim],Min[1:dim],Max[1:dim],N[1:dim],Len[1:dim],
    )
end

"""
    struct MeshCartesianStatic

## Fields
- `config::MeshConfig`
- `pos`
- `vel`
- `acc`
- `e`
- `rho`
- `phi`
- `B`
- `E`
- `U`
- `F`
- `G`
- `H`
- `J`

## Constructors
keywords are passed into `MeshConfig`
- `MeshCartesianStatic(config::MeshConfig, units = nothing)`
- `MeshCartesianStatic(units = nothing; kw...)`
- `MeshCartesianStatic(particles::StructArray, units = nothing; kw...)`
- `MeshCartesianStatic(particles::Array, units = nothing; kw...)`

For more instructions, see the documentation: https://juliaastrosim.github.io/PhysicalMeshes.jl/dev
"""
struct MeshCartesianStatic{I, VI, V, U, D, POS, VEL, ACC, _e, RHO, PHI, _B, _E, _U, _F, _G, _H, _J} <: AbstractMesh{U}
    config::MeshConfig{I,VI,V,U}
    data::D

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

function Base.show(io::IO, mesh::MeshCartesianStatic)
    n = fieldnames(typeof(mesh))
    print(io,
        """
        Static Cartesian Mesh

        $(mesh.config)

        Assigned field names: $(filter(x->!isnothing(getfield(mesh,x)), n[2:end]))
        """
    )
end

function __MeshCartesianStatic(config::MeshConfig, particles, ::VertexMode, units = nothing; gpu = false)
    a = [collect(config.Min[i] - config.Δ[i] * config.NG:config.Δ[i]:config.Max[i] + 1.000001*config.Δ[i] * config.NG) for i in 1:config.dim]

    iter = Iterators.product(a...)

    uLength = getuLength(units)

    zv = ZeroValue(units)

    if config.dim == 1 # StructArray is empty for eltype
        pos = [PVector(p..., uLength) for p in iter]
        vel = [zv.vel for p in iter]
        acc = [zv.acc for p in iter]
    else
        pos = StructArray(PVector(p..., uLength) for p in iter)
        vel = StructArray(zv.vel for p in iter)
        acc = StructArray(zv.acc for p in iter)
    end
    e = [zv.potpermass for p in iter]
    rho = [zv.density for p in iter]
    phi = [zv.potpermass for p in iter]

    if gpu
        return MeshCartesianStatic(
            config,
            cu(particles),
            cu(pos), cu(vel), cu(acc), cu(e), cu(rho), cu(phi),
            nothing, nothing,
            nothing, nothing, nothing, nothing, nothing,
        )
    else
        return MeshCartesianStatic(
            config,
            particles,
            pos, vel, acc, e, rho, phi,
            nothing, nothing,
            nothing, nothing, nothing, nothing, nothing,
        )
    end
end

function __MeshCartesianStatic(config::MeshConfig, particles, ::CellMode, units = nothing; kw...)

end

function MeshCartesianStatic(config::MeshConfig, particles, units = nothing; kw...)
    return __MeshCartesianStatic(config, particles, config.mode, units; kw...)
end

function MeshCartesianStatic(::Nothing, units = nothing; gpu = false, kw...)
    config = MeshConfig(units; kw...)
    return __MeshCartesianStatic(config, nothing, config.mode, units; gpu)
end

function MeshCartesianStatic(particles::StructArray, units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
    boundary = Periodic(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
    NG = 1,
    assign = true,
    gpu = false,
    kw...
)
    ratio = 1.01 + max(2/Nx, 2/Ny, 2/Nz) # Make sure that backward assignment are not affected by boundaries, at least in 3-element laplace conv
    e = extent(particles) * ratio

    config = MeshConfig(e, units; Nx, Ny, Nz, NG, mode, assignment, boundary, kw...)
    mesh = __MeshCartesianStatic(config, particles, mode, units; gpu)

    if assign
        assignmesh(particles, mesh, :Mass, :rho)
    end
    return mesh
end

function MeshCartesianStatic(particles::Array, units = nothing; kw...)
    return MeshCartesianStatic(StructArray(particles), units; kw...)
end