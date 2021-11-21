"""
    struct MeshConfig

- `mode`: the way of sampling physical properties
  - `CellMode`: properties are located in cell centers. `Nx × Ny × Nz` data points in total.
  - `VertexMode`: properties are located on grid points. `(Nx + 1) × (Ny + 1) × (Nz + 1)` data points in total.
- `Nx`, `Ny`, `Nz`: number of cells. In `VertexMode`, there are N+1 vertices in each direction
- `NG`: Number of boundary points
"""
struct MeshConfig{I,VI,V,U,F}
    mode::MeshMode
    assignment::MeshAssignment
    boundary::BoundaryCondition

    units::U
    dim::I
    NG::I # Number of ghost points. 

    # Vector info
    Δ::V
    Min::V
    Max::V
    N::VI    # length of data grid
    Len::VI  # total length (add two side of boundaries)
    LenH::VI # length with one side of boundary

    eps::F
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
                          eps: $(config.eps)
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
    eps = 1.0e-6,
)
    Δx = (xMax-xMin)/Nx
    Δy = (yMax-yMin)/Ny
    Δz = (zMax-zMin)/Nz

    Δ = SVector(Δx, Δy, Δz)
    Min = SVector(xMin, yMin, zMin)
    Max = SVector(xMax, yMax, zMax)
    N = SVector(Nx, Ny, Nz)
    Len = N .+ (2 * NG)
    LenH = N .+ NG
    return MeshConfig(
        mode,assignment,boundary,units,dim,NG,
        Δ[1:dim],Min[1:dim],Max[1:dim],N[1:dim],Len[1:dim],LenH[1:dim],
        eps
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
    eps = 1.0e-6,
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
    LenH = N .+ NG
    return MeshConfig(
        mode,assignment,boundary,units,dim,NG,
        Δ[1:dim],Min[1:dim],Max[1:dim],N[1:dim],Len[1:dim],LenH[1:dim],
        eps,
    )
end

struct MeshCartesianStatic{I, VI, V, U, F, POS, VEL, ACC, _e, RHO, PHI, _B, _E, _U, _F, _G, _H, _J} <: AbstractMesh{U}
    config::MeshConfig{I,VI,V,U,F}
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

function __MeshCartesianStatic(config::MeshConfig, ::VertexMode, units = nothing)
    a = [collect(config.Min[i] - config.Δ[i] * config.NG:config.Δ[i]:config.Max[i] + 1.000001*config.Δ[i] * config.NG) for i in 1:config.dim]

    iter = Iterators.product(a...)

    uLength = getuLength(units)

    zv = ZeroValue(units)

    if config.dim == 1 # StructArray is empty for eltype
        pos = SVector{config.N[1]+1+2*config.NG[1]}([PVector(p..., uLength) for p in iter])
        vel = SVector{config.N[1]+1+2*config.NG[1]}([zv.vel for p in iter])
        acc = SVector{config.N[1]+1+2*config.NG[1]}([zv.acc for p in iter])
    else
        pos = StructArray(PVector(p..., uLength) for p in iter)
        vel = StructArray(zv.vel for p in iter)
        acc = StructArray(zv.acc for p in iter)
    end
    e = MArray{Tuple{(config.N.+(1+2*config.NG))...}}([zv.potpermass for p in iter])
    rho = MArray{Tuple{(config.N.+(1+2*config.NG))...}}([zv.density for p in iter])
    phi = MArray{Tuple{(config.N.+(1+2*config.NG))...}}([zv.potpermass for p in iter])

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

function MeshCartesianStatic(units = nothing; kw...)
    config = MeshConfig(units; kw...)
    return __MeshCartesianStatic(config, config.mode, units)
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
    kw...
)
    ratio = 1.01 + max(2/Nx, 2/Ny, 2/Nz) # Make sure that backward assignment are not affected by boundaries, at least in 3-element laplace conv
    e = extent(particles) * ratio

    config = MeshConfig(e, units; Nx, Ny, Nz, NG, mode, assignment, boundary, kw...)
    mesh = __MeshCartesianStatic(config, mode, units)

    if assign
        assignmesh(particles, mesh, :Mass, :rho)
    end
    return mesh
end

function MeshCartesianStatic(particles::Array, units = nothing; kw...)
    return MeshCartesianStatic(StructArray(particles), units; kw...)
end