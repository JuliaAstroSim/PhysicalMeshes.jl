# Static Cartesian Mesh

```@repl example
using PhysicalMeshes, Unitful, UnitfulAstro # hide

# generate data
pos = [
    PVector(-1.0, -1.0, -1.0),
    PVector(-1.0, +1.0, -1.0),
    PVector(+1.0, -1.0, -1.0),
    PVector(+1.0, +1.0, -1.0),
    PVector(-1.0, -1.0, +1.0),
    PVector(-1.0, +1.0, +1.0),
    PVector(+1.0, -1.0, +1.0),
    PVector(+1.0, +1.0, +1.0),
];
dataArray = [Ball() for i in 1:8];
assign_particles(dataArray, :Pos, pos)
assign_particles(dataArray, :Mass, 1.0)

# construct mesh
m = MeshCartesianStatic(dataArray)

m.rho

# Check mass assignment
sum(m.rho) * m.config.Δ[1]^3 ≈ 8.0
```