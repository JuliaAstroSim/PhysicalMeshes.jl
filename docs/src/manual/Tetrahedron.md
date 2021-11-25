# Tetrahedron

```@repl example
using PhysicalMeshes, Unitful, UnitfulAstro # hide
t = Tetrahedron(
    PVector(1.0, 0.0, 1.0),
    PVector(1.0, 1.0, 0.0),
    PVector(0.0, 1.0, 1.0),
    PVector(1.0, 1.0, 1.0)
)
centroid(t)
circumcenter(t)
insphere(t, PVector(0.5, 0.5, 0.5))
insphere(t, PVector())
insphere(t, PVector(2.0, 0.0, 0.0))
PhysicalMeshes.volume(t)
orientation(t)
```