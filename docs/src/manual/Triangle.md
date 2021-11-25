# Triangle

## Triangle 2D

```@repl example
using PhysicalMeshes, Unitful, UnitfulAstro # hide

t = Triangle2D(PVector2D(4.0, 3.0), PVector2D(4.0, 0.0), PVector2D(0.0, 3.0))
orientation(t)
len(t)
PhysicalMeshes.area(t)
centroid(t)
circumcenter(t)
incircle(t, PVector2D(3.0, 2.0))
incircle(t, PVector2D())
incircle(t, PVector2D(6.0, 0.0))
```

## Triangle 3D

```@repl example
t = Triangle(PVector(0.0, 0.0, 0.0, u"m"), PVector(3.0, 4.0, 0.0, u"m"), PVector(3.0, 4.0, 12.0, u"m"))
len(t)
PhysicalMeshes.area(t)
centroid(t)
circumcenter(t)
```