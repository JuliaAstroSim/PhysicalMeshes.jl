using Base: IteratorEltype
abstract type MeshMode end
struct CellMode <: MeshMode end
struct VertexMode <: MeshMode end

"""
    struct MeshConfig

- `mode`: the way of sampling physical properties
  - `CellMode`: properties are located in cell centers. `Nx × Ny × Nz` data points in total.
  - `VertexMode`: properties are located on grid points. `(Nx + 1) × (Ny + 1) × (Nz + 1)` data points in total.
- `Nx`, `Ny`, `Nz`: number of cells 
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

function MeshConfig(;
        mode = VertexMode(),
        Nx = 10,
        Ny = 10,
        Nz = 10,
        xMin = -1.0,
        xMax = +1.0,
        yMin = -1.0,
        yMax = +1.0,
        zMin = -1.0,
        zMax = +1.0,
    )
        Δx = (xMax-xMin)/Nx
        Δy = (yMax-yMin)/Ny
        Δz = (zMax-zMin)/Nz
    return MeshConfig(mode,Nx,Ny,Nz,xMin,xMax,yMin,yMax,zMin,zMax,Δx,Δy,Δz)
end

struct MeshCartesianStatic{I, LEN, POS, VEL, ACC, E, RHO, PHI}
    config::MeshConfig{I, LEN}
    pos::POS
    vel::VEL
    acc::ACC
    e::E        # unit energy
    rho::RHO    # density
    phi::PHI    # potential

    U
    F
    G
    H
    J
end

function MeshCartesianStatic(;
        config = MeshConfig(),
    )
        x = collect(config.xMin:config.Δx:config.xMax)
        y = collect(config.yMin:config.Δy:config.yMax)
        z = collect(config.zMin:config.Δz:config.zMax)
        pos = PVector.(Iterators.product(x,y,z))
end