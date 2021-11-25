# Line

## 2D Line

```@repl example
using PhysicalMeshes, Unitful, UnitfulAstro # hide

line = Line2D(PVector2D(), PVector2D(1.0, 0.0))
len(line)
midpoint(line)
line + PVector(1.0,1.0)
line * 2
```

## 3D Space Line

```@repl example
line = Line(PVector(u"m"), PVector(1.0, 0.0, 0.0, u"m"))
len(line)
midpoint(line)
line + PVector(1.0,2.0,3.0,u"m")
line * 3
```