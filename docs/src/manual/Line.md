```@raw html
<!--
  Line manual page.
-->
```

# Line

`Line2D` and `Line` are 2-point segments, parameterised by their endpoints.
They support `+ / - / *scalar / /scalar / +vector / -vector` operations that
translate both endpoints simultaneously. The 2D version also has a unit
`normal(...)` (counter-clockwise convention).

## Types

The line type hierarchy is `AbstractGeometryType{T} > AbstractLine{T} >
AbstractLine{2D,3D}{T} > Line{2D}`. See [Library reference](@ref Methods) for
the concrete structs.

## Arithmetic

See [Library reference](@ref Methods) for `len`, `centroid`, `center`,
`midpoint`, `normal` on the line family.

## Example

```@repl line
using PhysicalMeshes
using PhysicalParticles, Unitful, UnitfulAstro

line_2d = Line2D(PVector2D(), PVector2D(1.0, 0.0))
len(line_2d)
midpoint(line_2d)
line_2d + PVector2D(1.0, 1.0)
line_2d * 2

line_3d = Line(PVector(u"m"), PVector(1.0, 0.0, 0.0, u"m"))
len(line_3d)
midpoint(line_3d)
line_3d + PVector(1.0, 2.0, 3.0, u"m")
line_3d * 3
```

## Properties and conventions

- `Line2D` <: `AbstractLine2D` <: `AbstractLine`. The constructor
  `Line(::AbstractPoint2D, ::AbstractPoint2D)` is provided as a convenience
  that auto-selects `Line2D` when the input is 2D.
- Arithmetic operators are explicitly noted: `Line + PVector` is
  *intentionally defined* (parallel translation of both endpoints), in
  contrast to the deliberate omission of `Line + Point` where the semantics
  would be ambiguous.
- `normal(::Line2D)` returns `normalize(perp(line.b - line.a))` — the
  counter-clockwise perpendicular, of the same unit as the endpoints.
- `centroid`, `center` and `midpoint` are aliases; they all return
  `(a + b) / 2`.
