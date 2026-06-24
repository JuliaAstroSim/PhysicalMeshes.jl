@testset "Triangle" begin
    @testset "Triangle2D" begin
        t = Triangle2D(PVector2D(4.0, 3.0), PVector2D(4.0, 0.0), PVector2D(0.0, 3.0))

        @test orientation(t) == NegativelyOriented()

        @test len(t) == 12.0
        @test PhysicalMeshes.area(t) == 6.0

        @test centroid(t) == PVector2D(8.0/3.0, 2.0)
        @test circumcenter(t) == PVector2D(2.0, 1.5)

        @test incircle(t, PVector2D(3.0, 2.0)) == Interior()
        @test incircle(t, PVector2D()) == OnEdge()
        @test incircle(t, PVector2D(6.0, 0.0)) == Exterior()
    end

    @testset "Triangle" begin
        t = Triangle(PVector(0.0, 0.0, 0.0, u"m"), PVector(3.0, 4.0, 0.0, u"m"), PVector(3.0, 4.0, 12.0, u"m"))

        @test len(t) == 30.0u"m"
        @test PhysicalMeshes.area(t) == 30.0u"m^2"

        @test centroid(t) == PVector(2.0, 8.0 / 3.0, 4.0, u"m")
        @test circumcenter(t) == PVector(1.5, 2.0, 6.0, u"m")
    end

    # 边界测试
    @testset "Degenerate Triangles" begin
        # 2D 退化三角形 - 三点共线
        t_degenerate_2d = Triangle2D(PVector2D(0.0, 0.0), PVector2D(1.0, 1.0), PVector2D(2.0, 2.0))
        @test PhysicalMeshes.area(t_degenerate_2d) == 0.0
        @test orientation(t_degenerate_2d) == UnOriented()

        # 3D 退化三角形 - 三点共线
        t_degenerate_3d = Triangle(PVector(0.0, 0.0, 0.0), PVector(1.0, 1.0, 1.0), PVector(2.0, 2.0, 2.0))
        @test PhysicalMeshes.area(t_degenerate_3d) == 0.0
    end

    @testset "Zero Area Triangles" begin
        # 两点重合
        t_zero_2d = Triangle2D(PVector2D(0.0, 0.0), PVector2D(0.0, 0.0), PVector2D(1.0, 0.0))
        @test PhysicalMeshes.area(t_zero_2d) == 0.0
        @test len(t_zero_2d) == 2.0  # 0 + 1 + 1

        # 三点重合
        t_all_same = Triangle2D(PVector2D(1.0, 1.0), PVector2D(1.0, 1.0), PVector2D(1.0, 1.0))
        @test PhysicalMeshes.area(t_all_same) == 0.0
        @test len(t_all_same) == 0.0
    end

    @testset "Orientation Tests" begin
        # 顺时针方向 (负方向)
        t_cw = Triangle2D(PVector2D(0.0, 0.0), PVector2D(0.0, 1.0), PVector2D(1.0, 0.0))
        @test orientation(t_cw) == NegativelyOriented()

        # 逆时针方向 (正方向)
        t_ccw = Triangle2D(PVector2D(0.0, 0.0), PVector2D(1.0, 0.0), PVector2D(0.0, 1.0))
        @test orientation(t_ccw) == PositivelyOriented()
    end

    @testset "Incircle Boundary Tests" begin
        # Triangle (0,0), (2,0), (0,2) has circumcenter (1, 1) and circumradius sqrt(2)
        t = Triangle2D(PVector2D(0.0, 0.0), PVector2D(2.0, 0.0), PVector2D(0.0, 2.0))

        # 三角形顶点恰好在 incircle 上 (距离 = sqrt(2) = 半径)
        @test incircle(t, PVector2D(0.0, 0.0)) == OnEdge()
        @test incircle(t, PVector2D(2.0, 0.0)) == OnEdge()
        @test incircle(t, PVector2D(0.0, 2.0)) == OnEdge()

        # (2,2) 也是 incircle 上的点
        @test incircle(t, PVector2D(2.0, 2.0)) == OnEdge()

        # 边上的点（不在 incircle 上） - 它们在 incircle 内部
        @test incircle(t, PVector2D(1.0, 0.0)) == Interior()
        @test incircle(t, PVector2D(0.0, 1.0)) == Interior()
        # (1,1) 是 incircle 圆心,在 incircle 内部
        @test incircle(t, PVector2D(1.0, 1.0)) == Interior()

        # 内部点
        @test incircle(t, PVector2D(0.5, 0.5)) == Interior()

        # 外部点
        @test incircle(t, PVector2D(3.0, 0.0)) == Exterior()
        @test incircle(t, PVector2D(0.0, 3.0)) == Exterior()
    end

    @testset "Large and Small Coordinates" begin
        # 极大坐标
        t_large = Triangle(PVector(1e10, 0.0, 0.0), PVector(0.0, 1e10, 0.0), PVector(0.0, 0.0, 1e10))
        @test PhysicalMeshes.area(t_large) > 0

        # 极小坐标
        t_small = Triangle(PVector(1e-10, 0.0, 0.0), PVector(0.0, 1e-10, 0.0), PVector(0.0, 0.0, 1e-10))
        @test PhysicalMeshes.area(t_small) > 0
    end

    @testset "Unitful Triangle Tests" begin
        # 带单位的三角形
        t_unit = Triangle(PVector(0.0, 0.0, 0.0, u"m"), PVector(3.0, 0.0, 0.0, u"m"), PVector(0.0, 4.0, 0.0, u"m"))
        @test PhysicalMeshes.area(t_unit) == 6.0u"m^2"
        @test len(t_unit) == 12.0u"m"

        # 验证 circumcenter 和 centroid 返回正确单位
        cc = circumcenter(t_unit)
        ce = centroid(t_unit)
        @test cc isa PVector
        @test ce isa PVector
        @test unit(cc.x) == u"m"
        @test unit(ce.x) == u"m"
    end
end
