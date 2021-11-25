# Cube

## 2D Cube

```@repl example
using PhysicalMeshes, Unitful, UnitfulAstro # hide

c = Cube(PVector2D(1.0, 1.0), PVector2D())
PhysicalMeshes.area(c)
```

## 3D Cube

```@repl example
c = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
PhysicalMeshes.volume(c)
```