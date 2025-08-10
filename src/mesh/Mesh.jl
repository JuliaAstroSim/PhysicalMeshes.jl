"""
$(TYPEDEF)
$(TYPEDFIELDS)

## Constructors
- `MeshConfig(units = nothing; kw...)`: all fields can be modified by keywords. Try construct one to see default values
- `MeshConfig(e::Extent, units = nothing; kw...)`: `Min` and `Max` are extracted from `e`. Other field can be modified by keywords

For more instructions, see the documentation: https://juliaastrosim.github.io/PhysicalMeshes.jl/dev
"""
struct MeshConfig{I,VI,V,U,D}
    """
    the way of sampling physical properties
    - `CellMode`: properties are located in cell centers. `Nx × Ny × Nz` data points in total.
    - `VertexMode`: properties are located on grid points. `(Nx + 1) × (Ny + 1) × (Nz + 1)` data points in total.
    """
    mode::MeshMode
    """
    mesh assignment algorithm
    - `NGP`: nearest grid point
    - `CIC`: cloud in cell
    - `TSC`: triangular shaped cloud
    """
    assignment::MeshAssignment
    """
    boundary conditions
    - `Periodic`
    - `Vacuum`
    - `Dirichlet`
    """
    boundary::BoundaryCondition
    "Enlarge the mesh"
    enlarge::Float64
    "Trait type of `DeviceType`. Supported: `CPU()`, `GPU()`. Defalut: `CPU()`"
    device::D

    "support `nothing`, `uAstro`, `uSI`, `uGadget2`, `uCGS`"
    units::U
    "dimension of the mesh"
    dim::I
    "Number of boundary ghost cells"
    NG::I

    # Vector info
    "mesh resolution in each direction"
    Δ::V
    "minimum coordinate of each direction, not including ghost cells. Corresponds to cell coordinates in cell mode, and vertex coordinates in vertex mode."
    Min::V
    "maximum coordinate of each direction, not including ghost cells. Corresponds to cell coordinates in cell mode, and vertex coordinates in vertex mode."
    Max::V
    "number of cells in each direction, not including ghost cells. In `VertexMode`, there are N+1 vertices in each direction"
    N::VI
    "total number of cells in each direction, including ghost cells."
    Len::VI
end

@inline length(p::T) where T <: MeshConfig = 1
@inline iterate(p::T) where T <: MeshConfig = (p,nothing)
@inline iterate(p::T,st) where T <: MeshConfig = nothing

function Base.show(io::IO, config::MeshConfig)
    print(io, 
    """
    Mesh config:
                          dim: $(config.dim)
                         mode: $(config.mode)
            assignment method: $(config.assignment)
           Boundary Condition: $(config.boundary)
           enlarge: $(config.enlarge)
           device: $(traitstring(config.device))
                        units: $(config.units)
           Number of Vertices: $(config.N.+1)
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
    yMin = xMin,
    yMax = xMax,
    zMin = xMin,
    zMax = xMax,
    dim = 3,
    device = CPU(),
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
        mode,assignment,boundary,
        1.0, device,
        units,dim,NG,
        Δ[1:dim],Min[1:dim],Max[1:dim],N[1:dim],Len[1:dim],
    )
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)

## Constructors
keywords are passed into `MeshConfig`
- `MeshCartesianStatic(config::MeshConfig, units = nothing)`
- `MeshCartesianStatic(units = nothing; kw...)`
- `MeshCartesianStatic(particles::StructArray, units = nothing; kw...)`
- `MeshCartesianStatic(particles::Array, units = nothing; kw...)`

For more instructions, see the documentation: https://juliaastrosim.github.io/PhysicalMeshes.jl/dev
"""
struct MeshCartesianStatic{I, VI, V, U, D, POS, VEL, ACC, _e, RHO, PHI, _B, _E, _j} <: AbstractMesh{U}
    config::MeshConfig{I,VI,V,U}
    data::D

    pos::POS
    vel::VEL
    acc::ACC
    "energy density"
    e::_e
    "density"
    rho::RHO
    "potential"
    phi::PHI

    # MHD
    "magnetic field"
    B::_B
    "eletric field"
    E::_E
    "charge density"
    rho_e
    "eletrical circuit field"
    j::_j
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

function __MeshCartesianStatic(config::MeshConfig, particles, ::VertexMode, units = nothing;
    mhd = false,
    data_on_cpu = false,
)
    a = [collect(LinRange(config.Min[i] - config.Δ[i] * config.NG, config.Max[i] + config.Δ[i] * config.NG, config.Len[i]+1)) for i in 1:config.dim]
    iter = Iterators.product(a...)

    zv = ZeroValue(eltype(ustrip(config.Δ[1])), units)

    if config.dim == 1 # StructArray is empty for eltype
        pos = [PVector(p...) for p in iter]
        vel = [zv.vel for p in iter]
        acc = [zv.acc for p in iter]
    else
        pos = StructArray(PVector(p...) for p in iter)
        vel = StructArray(zv.vel for p in iter)
        acc = StructArray(zv.acc for p in iter)
    end
    e = [zv.potpermass for p in iter]
    rho = [zv.density for p in iter]
    phi = [zv.potpermass for p in iter]

    if config.device isa GPU && !data_on_cpu
        return MeshCartesianStatic(
            config,
            cu(particles),
            cu(pos), cu(vel), cu(acc), cu(e), cu(rho), cu(phi),
            nothing, nothing, nothing, nothing,
        )
    else
        return MeshCartesianStatic(
            config,
            particles,
            pos, vel, acc, e, rho, phi,
            nothing, nothing, nothing, nothing,
        )
    end
end

function __MeshCartesianStatic(config::MeshConfig, particles, ::CellMode, units = nothing; kw...)

end

function MeshCartesianStatic(config::MeshConfig, particles, units = nothing; kw...)
    return __MeshCartesianStatic(config, particles, config.mode, units; kw...)
end

"""
$(TYPEDSIGNATURES)
Construct a static Cartesian mesh from nothing.

## Keywords
- `mhd::Bool`. If `true`, initiate `B`, `E`, `rho_e` and `j`. Default is `false`
"""
function MeshCartesianStatic(units::Union{Nothing, Vector{Unitful.FreeUnits{N, D, nothing} where {N, D}}} = nothing; mhd = false, data_on_cpu = false, kw...)
    config = MeshConfig(units; kw...)
    return __MeshCartesianStatic(config, nothing, config.mode, units; mhd, data_on_cpu)
end

"""
$(TYPEDSIGNATURES)
Construct a static Cartesian mesh containing particles.
The extent is enlarge by keyword argument `enlarge=2.01`.
"""
function MeshCartesianStatic(particles::StructArray, units = nothing;
    mode = VertexMode(),
    assignment = CIC(),
    boundary = Periodic(),
    Nx = 10,
    Ny = 10,
    Nz = 10,
    NG = 1,
    xMin = nothing,
    xMax = nothing,
    yMin = nothing,
    yMax = nothing,
    zMin = nothing,
    zMax = nothing,
    assign = true,
    device = CPU(),
    data_on_cpu = false,
    enlarge = 2.01,
    cube = true,
    kw...
)
    E = extent(particles) * enlarge
    if cube
        Min = min(E.xMin, E.yMin, E.zMin)
        !isnothing(xMin) && (Min = min(Min, xMin))
        !isnothing(yMin) && (Min = min(Min, yMin))
        !isnothing(zMin) && (Min = min(Min, zMin))
        Max = max(E.xMax, E.yMax, E.zMax)
        !isnothing(xMax) && (Max = max(Max, xMax))
        !isnothing(yMax) && (Max = max(Max, yMax))
        !isnothing(zMax) && (Max = max(Max, zMax))
        xMin = yMin = zMin = Min
        xMax = yMax = zMax = Max
    else
        xMin = isnothing(xMin) ? E.xMin : xMin
        xMax = isnothing(xMax) ? E.xMax : xMax
        yMin = isnothing(yMin) ? E.yMin : yMin
        yMax = isnothing(yMax) ? E.yMax : yMax
        zMin = isnothing(zMin) ? E.zMin : zMin
        zMax = isnothing(zMax) ? E.zMax : zMax
    end
    config = MeshConfig(units;
        Nx, Ny, Nz, NG,
        xMin, xMax, yMin, yMax, zMin, zMax,
        mode, assignment, boundary, device,
        kw...
    )
    mesh = __MeshCartesianStatic(config, particles, mode, units; data_on_cpu)

    if assign
        assignmesh(particles, mesh, :Mass, :rho)
    end
    return mesh
end

function MeshCartesianStatic(particles::Array, units = nothing; kw...)
    return MeshCartesianStatic(StructArray(particles), units; kw...)
end