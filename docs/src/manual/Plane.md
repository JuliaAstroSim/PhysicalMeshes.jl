```@raw html
<!--
  Plane manual page.
-->
```

# Plane

A `Plane{T}` is defined by three `PVector{T}` points (`a`, `b`, `c`). The
geometric structure is the triangle / disk that lies on the unique plane
through those points. Planes are used by [`Ray`](@ref) for ray-plane
intersection and by the [`Polygon`](@ref) constructors to check
coplanarity.

## Type

## Arithmetic and projection

See [Library reference](@ref Methods) for the `normal(::Point, ::Point, ::Point)`,
`normal(::Plane)`, and `distance(::AbstractPoint3D, ::Plane)` docstrings.

`distance(p, plane)` returns the **scalar** signed distance of the
point from the plane (without `abs`), so the sign encodes which side of
the plane the point lies on. Use `abs(...)` if you want the unsigned
distance.

## Coplanarity

See [Library reference](@ref Methods) for both `coplanar` methods
(point-vs-plane and vertex-collection vs plane).

`coplanar(p, plane, threshold)` returns `true` when the absolute
distance from `p` to `plane` is at most `threshold`. The base method
`coplanar(vertices, threshold)` is shared with [`Polygon3D`](@ref) and
checks that all of `vertices` lie on a single plane.

## Example

```@repl plane
using PhysicalMeshes
using PhysicalParticles

pl = Plane(PVector(0.0, 0.0, 0.0),
           PVector(1.0, 0.0, 0.0),
           PVector(0.0, 1.0, 0.0))

PhysicalMeshes.normal(pl)
PhysicalMeshes.distance(PVector(0.0, 0.0, 3.0), pl)
PhysicalMeshes.distance(PVector(0.0, 0.0, -3.0), pl)
PhysicalMeshes.coplanar(PVector(0.5, 0.5, 0.0), pl, 1e-9)
PhysicalMeshes.coplanar(PVector(0.5, 0.5, 1.0), pl, 1e-9)
```

## Properties and conventions

- `Plane` <: `AbstractPlane3D` <: `AbstractPlane`. There is currently no
  2D plane primitive; the equivalent construct in 2D is a `Line2D`.
- The normal returned by `normal(...)` has the **same length units** as the
  input points — i.e. for unitful input it is `m`, not `m/m`. This means
  that downstream dot products preserve the input's length² units (e.g.
  `m^2` for area computations).
