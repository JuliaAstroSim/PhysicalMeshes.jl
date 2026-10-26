```@raw html
<!--
  Sphere manual page.
-->
```

# Sphere

`Sphere` is a minimal 3D primitive: a center point plus an (unsigned)
radius. Unit-aware (`radius` carries the same units as the center). Used
externally by the JuliaAstroSim ecosystem for gravitational softening
checks and surface-boundary tests.

## Type

## Membership

See [Library reference](@ref Methods) for `interior(::Sphere, ::Any)` and
`exterior(::Sphere, ::Any)`.

## Example

```@repl sphere
using PhysicalMeshes
using PhysicalParticles, Unitful, UnitfulAstro

s = Sphere(PVector(0.0, 0.0, 0.0, u"kpc"), 1.0u"kpc")
interior(s, PVector(0.5, 0.0, 0.0, u"kpc"))    # true
interior(s, PVector(1.5, 0.0, 0.0, u"kpc"))    # false
exterior(s, PVector(2.0, 0.0, 0.0, u"kpc"))    # true
exterior(s, PVector(0.5, 0.0, 0.0, u"kpc"))    # false
```

## Properties and conventions

- `interior` uses strict (`<`), `exterior` uses strict (`>`). A point
  exactly on the sphere surface returns `false` for both; use
  `isapprox(norm(p - s.center), s.radius; atol=...)` for the closed-interval
  version.
- `Sphere.center` and `Sphere.radius` are the only fields; radius is a
  scalar with the same units as the center.
