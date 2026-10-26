```@raw html
<!--
  Guide / landing page for the Manual section.
-->
```

# Guide

```@contents
Pages = ["guide.md"]
Depth = 3
```

`PhysicalMeshes.jl` covers three layers of functionality. The first three
sections of this guide walk through each one in turn with runnable examples.
The later sections cover the parallel / halo-exchange code path and the
solver utilities that round out the package.

## Installation

From the Julia REPL, enter the Pkg REPL with `]`:

```julia
pkg> add PhysicalMeshes
```

Or install from GitHub to track `master`:

```julia
pkg> add https://github.com/JuliaAstroSim/PhysicalMeshes.jl
```

To verify the install:

```julia
julia> using PhysicalMeshes

julia> Triangle(PVector(0.0, 0.0), PVector(1.0, 0.0), PVector(0.0, 1.0)) |> PhysicalMeshes.area
0.5
```

To run the upstream test suite (covers geometry predicates, mesh construction
and halo exchange):

```julia
pkg> test PhysicalMeshes
```

## Coordinate types and units

`PhysicalMeshes.jl` is built on `PhysicalParticles.jl`, which provides the
two point types used everywhere:

- `PVector2D{T}` — a 2-component `T`. Use it for 2D problems.
- `PVector{T}` — a 3-component `T`. Use it for 3D problems.

`T` can be either `Float64` (dimensionless) or `Unitful.Quantity`. Every
arithmetic predicate and constructor preserves units.

```@example units
using PhysicalMeshes
using PhysicalParticles
using Unitful, UnitfulAstro

# Plain (dimensionless) points
dim_less = Triangle2D(PVector2D(0.0, 0.0),
                      PVector2D(1.0, 0.0),
                      PVector2D(0.0, 1.0))
PhysicalMeshes.area(dim_less)            # 0.5

# Unitful points
with_units = Triangle(PVector(0.0, 0.0, 0.0, u"kpc"),
                      PVector(3.0, 0.0, 0.0, u"kpc"),
                      PVector(0.0, 4.0, 0.0, u"kpc"))
PhysicalMeshes.area(with_units)         # 6 kpc^2
```

`uAstro` is a `Unitful` unit module bundled with `UnitfulAstro.jl` and is
the recommended unit list when constructing `MeshConfig` for astrophysics
runs (`uSI`, `uCGS`, `uGadget2` are also accepted).

## Geometry primitives

The geometry layer is intentionally small. Each primitive has:

- `len(...)` — perimeter (or length, for a line) — preserves units.
- `area(...)` / `PhysicalMeshes.volume(...)` — area / volume.
- `centroid(...)`, `midpoint(...)`, `normal(...)` — basic vector quantities.
- `orient(...)`, `orientation(...)` — signed orientation.
- `circumcenter(...)` — circumcenter / circumsphere center.
- `incircle(...)` / `insphere(...)` — robust Shewchuk-style interior /
  exterior / on-edge tests.
- `interior(...)` / `exterior(...)` — half-open point-membership predicates.

See the per-primitive pages for the full set:

- [`Line`](@ref), [`Plane`](@ref), [`Triangle`](@ref), [`Cube`](@ref),
  [`Tetrahedron`](@ref)
- [`Polygon`](@ref), [`Sphere`](@ref), [`Ray`](@ref)

### Robust predicates

`incircle(::Triangle2D, ::PVector2D)` and
`insphere(::Tetrahedron, ::PVector3D)` combine the orientation kernel and the
in-circle / in-sphere kernel so that the answer is invariant to whether the
input triangle is CCW or CW. The return type is one of the trait singletons
`Interior()`, `Exterior()`, or `OnEdge()`.

```@example predicate
using PhysicalMeshes
using PhysicalParticles

t = Triangle2D(PVector2D(4.0, 3.0),
               PVector2D(4.0, 0.0),
               PVector2D(0.0, 3.0))

incircle(t, PVector2D(3.0, 2.0))    # Interior()
incircle(t, PVector2D())            # OnEdge()
incircle(t, PVector2D(6.0, 0.0))    # Exterior()
```

## Building meshes

```@example mesh
using PhysicalMeshes
using PhysicalParticles
using AstroSimBase: Periodic, CPU
using Unitful, UnitfulAstro

# Bare config — no particles attached
cfg = MeshConfig(uAstro; Nx = 5, Ny = 5, Nz = 5,
                 mode  = VertexMode(),
                 assignment = CIC(),
                 boundary   = Periodic(),
                 device     = CPU())

m = MeshCartesianStatic(cfg)
m.rho                   # ArrayScalarField{Float64} of size (6, 6, 6)
```

The size of `m.rho` is `(Nx+1, Ny+1, Nz+1)` in `VertexMode` because there are
`N + 1` vertices along each axis. Switch to `CellMode` to get a `(Nx, Ny,
Nz)` field instead.

The full constructor accepts a wide range of keyword arguments. The most
relevant are:

| Keyword     | Meaning                                                       |
|-------------|---------------------------------------------------------------|
| `mode`      | `VertexMode()` (default) or `CellMode()`                      |
| `assignment`| `CIC()` (default) / `NGP()` / `TSC()` mass-assignment kernel |
| `boundary`  | `Periodic()` (default) / `Dirichlet()` / `Vacuum()`           |
| `device`    | `CPU()` (default) / `GPU()`                                   |
| `mhd`       | if `true`, allocate `B`, `E`, `rho_e`, `j` as non-zero buffers |
| `NG`        | ghost-cell width in each direction (default `1` when constructing from particles, `0` otherwise) |
| `enlarge`   | factor by which to enlarge the auto-derived extent (default `2.01`) |
| `assign`    | whether to call `assignmesh(particles, m, :Mass, :rho)` in the constructor (default `true`) |
| `units`     | positional; one of `uAstro`, `uSI`, `uCGS`, `uGadget2` or `nothing` |

See the full list in [`MeshConfig`](@ref).

## Assignment schemes

`PhysicalMeshes.jl` ships with three mass-assignment kernels, all of which
work on a `MeshCartesianStatic` whose `config.assignment` matches:

| Kernel    | Stencil size              | Properties                                                        |
|-----------|---------------------------|-------------------------------------------------------------------|
| `NGP()`   | 1 cell                    | Cheapest, no smoothing, mass not exactly conserved               |
| `CIC()`   | 2 cells per active axis   | (default) Mass-conserving, smooth, cheap                         |
| `TSC()`   | 3 cells per active axis   | Smoother than CIC, slightly more expensive                       |

To switch assignment mode:

```julia
m = MeshCartesianStatic(pos, uAstro;
                        Nx = 4, Ny = 4, Nz = 4,
                        assignment = TSC())
```

`assignmesh(particles, mesh, :Mass, :rho)` scatters the `:Mass` field of
`particles` onto `mesh.rho`; `assignparticle(particles, mesh, :Acc, :acc)`
interpolates in the reverse direction from `mesh.acc` onto a `:Acc` field
of `particles`.

Both calls dispatch on `mesh.config.assignment` and `mesh.config.mode`; you
do not have to specify them explicitly.

## Fields

Each mesh owns:

| Name    | Type                  | Meaning                                          |
|---------|-----------------------|--------------------------------------------------|
| `pos`   | `StructArray{PVector}` (≥2D) / `Vector{PVector}` (1D) | Vertex positions                  |
| `vel`   | vector field          | Initialized to zero                              |
| `acc`   | vector field          | Initialized to zero                              |
| `e`     | `ArrayScalarField`    | Energy density                                   |
| `rho`   | `ArrayScalarField`    | Mass density                                     |
| `phi`   | `ArrayScalarField`    | Gravitational / electrostatic potential          |
| `B`     | `ArrayVectorField`    | Magnetic field (MHD)                             |
| `E`     | `ArrayVectorField`    | Electric field (MHD)                             |
| `rho_e` | `ArrayScalarField`    | Charge density (PIC)                             |
| `j`     | `ArrayVectorField`    | Current density (PIC)                            |

The fields are wrapped in `ArrayScalarField` /
[`ArrayVectorField`](@ref) / [`ArrayTensorField`](@ref) types so that the
mass-assignment and halo-exchange code can dispatch on a unified interface
without having to know whether the underlying storage is `Array`,
`CuArray`, or `DArray`. Common reductions (`sum`, `prod`, `minimum`,
`maximum`, `extrema`) and broadcasts work directly on the wrappers. For
slice-assignment via `m.phi[2:end-1, ...] .= rhs`, the wrappers forward the
broadcast destination through `Base.Broadcast.dotview` so the assignment
actually updates the underlying data.

To get the raw storage, use `.data` (the field itself behaves like the
underlying `Array` in indexing and arithmetic).

## Parallel backends

`PhysicalMeshes.jl` ships with a backend abstraction that hides whether your
work lives on the host, on a `CuArray`, or on a `DArray`. The two entry
points are `detect_resources` and `select_backend`:

```@example backend
using PhysicalMeshes

resources = detect_resources()                # NamedTuple of counts
backend   = select_backend()                 # picks the best available
```

The returned `backend::ParallelBackend` carries the per-kind closures
(`device_array`, `distribute_array`, `local_array`, `release!`) and the
selected `fft_threads` count. Pass it to the wrapper functions listed in
the [`Parallel`](@ref manual/Parallel) page instead of branching on
`backend.kind` in user code.

## Distribution and halo exchange

`DistributedMesh` partitions a global `MeshConfig` across worker processes
along the `x`-axis (1D), a `rows × cols` grid (2D), or a slab along `x` (3D).
After construction, every worker owns a local `MeshCartesianStatic` plus
the list of neighbours needed for halo exchange.

A typical iteration looks like:

```julia
# On every worker
mesh = DistributedMesh(cfg, uAstro)

# ... push particles, deposit charge, solve fields, etc. ...

# Synchronise E, B, rho_e, etc. across neighbours
sync_mesh(mesh)
```

See the [`Distributed Mesh`](@ref manual/Distributed-Mesh) page for the
full protocol and the underlying `exchange_boundary_data` rendezvous.

## Algorithms

`PhysicalMeshes.jl` includes two reference solvers that demonstrate how to
build on top of the mesh primitives:

- `FDTD{1,2,3}D` — finite-difference time-domain electromagnetics.
- `PIC{1,2,3}D` — particle-in-cell, with `boris_algorithm` for velocity
  pushing, CIC-based `deposit_charge` / `deposit_current`, and a Jacobi
  `relax_poisson` field solver.

These steppers are meant as **references** — production runs typically use
the more specialized implementations in `PhysicalPlasma.jl` /
`NonlinearOptics.jl`. They do, however, exercise every field, every
boundary condition, and the halo exchange code path end-to-end.

```julia
fdtd = FDTD3D()
fdtd_simulate(m, fdtd, 1e-4, 100)        # 100 FDTD steps, dt = 1e-4

pic = PIC3D()
pic_simulate(m, pic, 1e-3, 50)           # 50 PIC steps
```

See the [`Algorithms`](@ref manual/Algorithms) page for details.

## Mesh tools

`zoom(a, ratio)` and `shrink(a, ratio)` resize 1D / 2D / 3D arrays by an
integer ratio. `zoom` accepts an `Interpolations` interpolator as a keyword
to do smooth interpolation instead of nearest-neighbour replication:

```@example resize
using PhysicalMeshes
using Interpolations

zoom([1, 5, 3], 2, interp = LinearInterpolation)
# 6-element Vector{Float64}:
# 1.0, 2.6, 4.2, 4.6, 3.8, 3.0

shrink(ones(6, 6, 6), 3)
# 2×2×2 Array{Float64, 3}: filled with 1.0
```

For mesh-specific reductions:

- `mass_center(m::MeshCartesianStatic)` returns the density-weighted mean
  position, falling back to the mean vertex position when `rho` is empty.
- `total_mass(m::MeshCartesianStatic)` integrates `m.rho` over the cell
  volume `prod(config.Δ)`.

## Visualisation

`PhysicalMeshes.jl` does not bundle its own plotting routines. The
companion package [`AstroPlot.jl`](https://github.com/JuliaAstroSim/AstroPlot.jl)
provides Makie / GLMakie / CairoMakie recipes for all the primitives and
mesh types documented here.
