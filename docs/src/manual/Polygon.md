```@raw html
<!--
  Polygon manual page.
-->
```

# Polygon

`Polygon2D` and `Polygon3D` carry an `Vector{PVector{...}}` of vertices
(>2 vertices, otherwise the constructor throws). Both support robust
point-in-polygon checks via ray-casting, convex-polygon tests, and signed
area / coplanarity tests.

## Types

The polygon type hierarchy is `AbstractGeometryType{T} > AbstractPolygon{T} >
AbstractPolygon{2D,3D}{T} > Polygon{2D,3D}`. See [Library reference](@ref Methods)
for the concrete structs.

## Factory helpers

See [Library reference](@ref Methods) for `polygon_rect`, `polygon_regular`.

## Area, coplanarity, convexity

See [Library reference](@ref Methods) for `area`, `normal`, `isconvex`,
`coplanar` and the related shape predicates.

## Point membership

See [Library reference](@ref Methods) for the two `is_inbound(::PVector, ::Polygon)`
methods.

For 2D inputs, `is_inbound` picks the convex (cross-product) algorithm
when the polygon is convex and falls back to the standard ray-casting
algorithm otherwise.

For 3D inputs, the point is first checked against the plane of the
polygon (with a small `unit(p.x)^2 * 1e-6` threshold) and then the 2D test
is run on a local coordinate system constructed from the first two edges.

## Compatibility shim

For backward compatibility with downstream code, both `Polygon2D` and
`Polygon3D` expose a `points` alias for their canonical `vertices` field
(through `Base.getproperty` and `Base.propertynames`).

## Example

```@repl polygon
using PhysicalMeshes
using PhysicalParticles, Unitful

# Rectangle
rect = polygon_rect(4.0, 2.0)
PhysicalMeshes.area(rect)        # 8.0
isconvex(rect)                   # true
PhysicalMeshes.is_inbound(PVector2D(0.0, 0.0), rect)   # true
PhysicalMeshes.is_inbound(PVector2D(5.0, 0.0), rect)   # false

# Regular hex (circumradius = 1.0)
hex = polygon_regular(6, 1.0)
PhysicalMeshes.area(hex)         # ≈ 2.598
isconvex(hex)                    # true
PhysicalMeshes.is_inbound(PVector2D(0.0, 0.0), hex)   # true
PhysicalMeshes.is_inbound(PVector2D(2.0, 0.0), hex)   # false
```

## Properties and conventions

- The internal canonical field is `vertices`. Code written against the
  older `.points` alias still works through a `Base.getproperty` shim.
- `area(::Polygon2D)` returns an unsigned value via the shoelace formula
  (decomposed as a fan of triangles anchored at the first vertex for
  numerical robustness).
- `normal(::Polygon2D)` returns the **unitless** `ẑ` direction (and
  preserves the input units); `normal(::Polygon3D)` calls the
  `normal(a, b, c)` three-point version.
