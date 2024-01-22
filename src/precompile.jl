@setup_workload begin
    @compile_workload begin
        # Tools
        zoom([1],2)
        zoom(ones(Int, 1, 1),2)
        zoom(ones(Int, 1, 1, 1),2)

        a1 = [1,2,3,4]
        a2 = a1*a1'
        a3 = ones(Int,4,4,4) .* reshape(a1,1,1,4)
        shrink(a1,2)
        shrink(a2,2)
        shrink(a3,2)

        # Basic types
        c = Cube(PVector2D(1.0, 1.0), PVector2D())
        PhysicalMeshes.area(c)
        c = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
        PhysicalMeshes.volume(c)
        interior(c, PVector(0.5, 0.5, 0.5, u"m"))
        exterior(c, PVector(1.5, 1.5, 1.5, u"m"))

        line = Line2D(PVector2D(), PVector2D(1.0, 0.0))
        len(line)
        midpoint(line)
        m = Line2D(PVector2D(2.0, 2.0, u"m"), PVector2D(2.0, 2.0, u"m"))
        n = Line2D(PVector2D(1.0, 1.0, u"m"), PVector2D(1.0, 1.0, u"m"))
        m + n
        m - n
        n * 2
        2 * n
        n + PVector2D(1.0, 1.0, u"m")
        m - PVector2D(1.0, 1.0, u"m")
        m / 2

        line = Line(PVector(u"m"), PVector(1.0, 0.0, 0.0, u"m"))
        len(line)
        midpoint(line)
        m = Line(PVector(2.0, 2.0, 2.0), PVector(2.0, 2.0, 2.0))
        n = Line(PVector(1.0, 1.0, 1.0), PVector(1.0, 1.0, 1.0))
        m + n
        m - n
        n * 2
        2 * n
        n + PVector(1.0, 1.0, 1.0)
        m - PVector(1.0, 1.0, 1.0)
        m / 2

        s = Sphere(PVector(), 1.0)
        interior(s, PVector(0.1, 0.1, 0.1))
        exterior(s, PVector(1.0, 1.0, 1.0))

        t = Tetrahedron(
            PVector(1.0, 0.0, 1.0),
            PVector(1.0, 1.0, 0.0),
            PVector(0.0, 1.0, 1.0),
            PVector(1.0, 1.0, 1.0)
        )
        centroid(t)
        circumcenter(t)
        insphere(t, PVector(0.5, 0.5, 0.5))
        insphere(t, PVector())
        insphere(t, PVector(2.0, 0.0, 0.0))
        PhysicalMeshes.volume(t)
        orientation(t)

        t = Triangle2D(PVector2D(4.0, 3.0), PVector2D(4.0, 0.0), PVector2D(0.0, 3.0))
        orientation(t)
        len(t)
        PhysicalMeshes.area(t)
        centroid(t)
        circumcenter(t)
        incircle(t, PVector2D(3.0, 2.0))
        incircle(t, PVector2D())
        incircle(t, PVector2D(6.0, 0.0))

        t = Triangle(PVector(0.0, 0.0, 0.0, u"m"), PVector(3.0, 4.0, 0.0, u"m"), PVector(3.0, 4.0, 12.0, u"m"))
        len(t)
        PhysicalMeshes.area(t)
        centroid(t)
        circumcenter(t)


        # Mesh
        pos = [
            PVector(-1.0, -1.0, -1.0),
            PVector(-1.0, +1.0, -1.0),
            PVector(+1.0, -1.0, -1.0),
            PVector(+1.0, +1.0, -1.0),
            PVector(-1.0, -1.0, +1.0),
            PVector(-1.0, +1.0, +1.0),
            PVector(+1.0, -1.0, +1.0),
            PVector(+1.0, +1.0, +1.0),
        ]
        dataArray = [Ball() for i in 1:8]
        assign_particles(dataArray, :Pos, pos)
        assign_particles(dataArray, :Mass, 1.0)
        m = MeshCartesianStatic(dataArray)

        #TODO Unitful mesh
    end
end
