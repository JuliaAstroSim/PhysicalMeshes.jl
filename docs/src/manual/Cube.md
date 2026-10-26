```@raw html
<!--
  Cube manual page.
-->
```

# Cube

`Cube2D` and `Cube` are axis-aligned bounding boxes defined by `top` and
`below` corner points. They are used internally by `MeshCartesianStatic`
when constructing a mesh from an `extent(particles)` and externally by
the rest of the JuliaAstroSim ecosystem for early-rejection collision
checks and AABB hierarchies.

## Types

The cube type hierarchy is `AbstractGeometryType{T} > AbstractCube{T} >
AbstractCube{2D,3D}{T} > Cube{2D}`. See [Library reference](@ref Methods) for
the concrete structs.

## Area, volume, membership

See [Library reference](@ref Methods) for `area`, `volume`, `interior`,
`exterior` on the cube family.

## Example

```@repl cube
using PhysicalMeshes
using PhysicalParticles, Unitful, UnitfulAstro

c_2d = Cube2D(PVector2D(1.0, 1.0), PVector2D())
PhysicalMeshes.area(c_2d)

c_3d = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
PhysicalMeshes.volume(c_3d)
```

## Properties and conventions

- The constructor `Cube(top::PVector2D, below::PVector2D)` is provided as a
  convenience that selects `Cube2D` automatically.
- `interior` and `exterior` use **strict** inequalities: points exactly on
  the cube's edge are reported as `false` for `interior` and `false` for
  `exterior`. If you need closed-interval semantics, test against
  `area == 0` (degenerate cases) or use `MeshConfig.Min / Max` directly.
- `area` and `volume` are **signed** quantities:
  `area(cube) = (top.x - below.x) * (top.y - below.y)`. Callers that need
  a non-negative result should wrap it in `abs(...)` — the package does not
  do this automatically to preserve orientation parity for downstream code.
