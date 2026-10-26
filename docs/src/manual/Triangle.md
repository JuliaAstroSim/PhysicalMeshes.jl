```@raw html
<!--
  Triangle manual page.
-->
```

# Triangle

`Triangle` (3D, generic) and `Triangle2D` are the work-horses of the
geometry layer: a triangle is enough topology to test orientation, circumcenter
and incircle membership, and the 2D / 3D code paths share the same algorithmic
families (Shewchuk-style robust predicates).

## Types

The triangle type hierarchy is `AbstractGeometryType{T} > AbstractTriangle{T}
> AbstractTriangle{2D,3D}{T} > Triangle{2D}`. See [Library reference](@ref Methods)
for the concrete structs.

## Arithmetic and basic vector quantities

See [Library reference](@ref Methods) for the canonical docstrings of
`len`, `orient`, `orientation`, `area`, `centroid`, `circumcenter`,
`incircle` on the triangle family.

## Example

```@repl triangle
using PhysicalMeshes
using PhysicalParticles, Unitful, UnitfulAstro

t2 = Triangle2D(PVector2D(4.0, 3.0),
                PVector2D(4.0, 0.0),
                PVector2D(0.0, 3.0))

orientation(t2)
len(t2)
PhysicalMeshes.area(t2)
centroid(t2)
circumcenter(t2)
incircle(t2, PVector2D(3.0, 2.0))
incircle(t2, PVector2D())
incircle(t2, PVector2D(6.0, 0.0))

t3 = Triangle(PVector(0.0, 0.0, 0.0, u"m"),
              PVector(3.0, 4.0, 0.0, u"m"),
              PVector(3.0, 4.0, 12.0, u"m"))

len(t3)
PhysicalMeshes.area(t3)
centroid(t3)
circumcenter(t3)
```

## Properties and conventions

- The 3D `Triangle` carries the units of its vertices — `len` returns a
  length, `area` an area², `centroid` and `circumcenter` a `PVector` with the
  same length units.
- `orientation(::Triangle2D)` returns one of the singletons
  `PositivelyOriented()`, `NegativelyOriented()` or `UnOriented()`. The same
  convention is used for [`Tetrahedron`](@ref) and other orientation-aware
  primitives.
- `incircle(::Triangle2D, p)` reports whether `p` lies inside, on, or
  outside the triangle's circumcircle. The test is Shewchuk-style and is
  invariant under reflection of the input triangle.
