```@raw html
<!--
  Ray manual page.
-->
```

# Ray

`Ray2D` and `Ray3D` represent a parametric ray `x + t * n`, `t ≥ 0`. Ray
geometry powers the laser propagation examples in the JuliaAstroSim
ecosystem and the [Finite-Difference Time-Domain](@ref manual/Algorithms)
electromagnetics solver. Each ray carries an `origin` (`x`) and a
`direction` (`n`); the 2D version additionally stores the azimuth angle
`theta` for cheap rotation / coordinate conversions.

## Types

The ray type hierarchy is `AbstractGeometryType{T} > AbstractRay{T} >
AbstractRay{2D,3D}{T} > Ray{2D}`. See [Library reference](@ref Methods) for
the concrete structs.

## Intersections and reflections

See [Library reference](@ref Methods) for `intersect` / `reflect` on the
ray family (and the generic [`reflect(v, n)`](@ref) vector mirror).

The package-level `reflect(vec, n)` mirrors a vector `vec` about a unit
normal `n`: `vec − 2·(n · vec) · n`.

## Example

```@repl ray
using PhysicalMeshes
using PhysicalParticles

ray = Ray2D(PVector2D(0.0, 0.5), PVector2D(1.0, 0.0))
line = Line2D(PVector2D(0.5, 0.0), PVector2D(0.5, 1.0))

hit, point = intersect(ray, line)
hit
point
reflected = reflect(ray, line)
reflected

# 3D ray-plane
r3 = Ray3D(PVector(0.0, 0.0, 0.0), PVector(1.0, 1.0, -1.0))
pl = Plane(PVector(1.0, 0.0, 0.0), PVector(0.0, 1.0, 0.0), PVector(0.0, 0.0, 1.0))
hit3, p3 = intersect(r3, pl)
hit3
p3
```

## Properties and conventions

- `Ray2D` zero / non-zero checks are dispatching on the **direction
  vector**: `zero(::Ray2D)` produces a degenerate ray with both origin and
  direction zeroed and `iszero(::Ray2D)` is true exactly when `iszero(r.n)`.
  The same convention holds for `Ray3D`.
- `intersect(::Ray2D, ::Line2D)` returns `(false, PVector2D())` for
  parallel rays and degenerate segments (`a == b`); otherwise it returns
  `(true, intersection_point)`.
- `intersect(::Ray3D, ::Plane)` treats the plane normal as the bisector:
  when `dot(ray.n, N) ≈ 0` the ray is parallel to the plane and we return
  `(false, PVector(0, 0, 0))`. A parallel-but-ill-defined ray is signalled
  by the `1e-6 * unit(denom)` tolerance rather than exact zero.
- The azimuth angle `theta` stored in `Ray2D` is automatically computed
  from the direction vector at construction; mutating `r.n` does **not**
  update it.
