```@raw html
<!--
  Tetrahedron manual page.
-->
```

# Tetrahedron

`Tetrahedron{T}` is the simplest 3D simplex: four vertices and the
corresponding circumcenter / insphere / volume / orientation machinery.

## Types

The tetrahedron type hierarchy is `AbstractGeometryType{T} >
AbstractTetrahedron{T} > Tetrahedron{T}`. See [Library reference](@ref Methods)
for the concrete struct.

## Arithmetic, orientation, circumcenter, insphere

See [Library reference](@ref Methods) for `centroid`, `circumcenter`,
`insphere`, `orient`, `volume`, `orientation` on the tetrahedron family.

## Example

```@repl tetra
using PhysicalMeshes
using PhysicalParticles

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

## Properties and conventions

- `volume(::Tetrahedron)` returns `abs(orient) / 6`, where `orient` is the
  signed 6×volume of the tetrahedron.
- `insphere(...)` uses the Shewchuk kernel and rejects near-coplanar
  configurations by combining the kernel value with the (orientation-aware)
  determinant of the four vertices. The return value is one of
  `Interior()`, `Exterior()`, or `OnEdge()`.
- `circumcenter(::Tetrahedron)` throws `ArgumentError("Coplanar points:
  no unique circumcenter")` for the degenerate 4-coplanar-points case.
