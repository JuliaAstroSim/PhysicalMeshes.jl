@testset "Tetrahedron" begin
    t = Tetrahedron(
        PVector(1.0, 0.0, 1.0),
        PVector(1.0, 1.0, 0.0),
        PVector(0.0, 1.0, 1.0),
        PVector(1.0, 1.0, 1.0)
    )

    @test centroid(t) == PVector(3.0, 3.0, 3.0) / 4.0
    @test circumcenter(t) == PVector(0.5, 0.5, 0.5)
    
    @test insphere(t, PVector(0.5, 0.5, 0.5)) == Interior()
    @test insphere(t, PVector()) == OnEdge()
    @test insphere(t, PVector(2.0, 0.0, 0.0)) == Exterior()

    @test PhysicalMeshes.volume(t) == 1.0 / 6.0

    @test orientation(t) == PositivelyOriented()
end

@testset "Tetrahedron - Degenerate (coplanar points)" begin
    # 四点共面 - 退化四面体（所有点 z=0）
    t_degenerate = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(1.0, 0.0, 0.0),
        PVector(0.0, 1.0, 0.0),
        PVector(1.0, 1.0, 0.0)
    )
    @test PhysicalMeshes.volume(t_degenerate) == 0.0
    @test orientation(t_degenerate) == UnOriented()
    
    # 另一种共面情况（点在 x-y 平面）
    t_degenerate2 = Tetrahedron(
        PVector(0.0, 0.0, 1.0),
        PVector(1.0, 0.0, 1.0),
        PVector(0.0, 1.0, 1.0),
        PVector(1.0, 1.0, 1.0)
    )
    @test PhysicalMeshes.volume(t_degenerate2) == 0.0
    @test orientation(t_degenerate2) == UnOriented()
end

@testset "Tetrahedron - Zero volume" begin
    # 三点重合
    t_zero_1 = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(0.0, 0.0, 0.0),
        PVector(0.0, 0.0, 0.0),
        PVector(1.0, 0.0, 0.0)
    )
    @test PhysicalMeshes.volume(t_zero_1) == 0.0
    @test orientation(t_zero_1) == UnOriented()
    @test centroid(t_zero_1) == PVector(0.25, 0.0, 0.0)

    # 四点重合
    t_zero_2 = Tetrahedron(
        PVector(1.0, 1.0, 1.0),
        PVector(1.0, 1.0, 1.0),
        PVector(1.0, 1.0, 1.0),
        PVector(1.0, 1.0, 1.0)
    )
    @test PhysicalMeshes.volume(t_zero_2) == 0.0
    @test centroid(t_zero_2) == PVector(1.0, 1.0, 1.0)
    @test orientation(t_zero_2) == UnOriented()
    
    # 两点重合
    t_zero_3 = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(1.0, 0.0, 0.0),
        PVector(0.0, 1.0, 0.0),
        PVector(0.0, 0.0, 0.0)  # 与第一点重合
    )
    @test PhysicalMeshes.volume(t_zero_3) == 0.0
    @test orientation(t_zero_3) == UnOriented()
end

@testset "Tetrahedron - Insphere boundary" begin
    # 4 个顶点均在 insphere 上 (OnEdge); 面上非顶点距离在 circumsphere 内部;
    # (1,1,1) 距 circumsphere 中心 (0.5,0.5,0.5) 为 sqrt(0.75) = 半径,恰好在表面上
    t = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(1.0, 0.0, 0.0),
        PVector(0.0, 1.0, 0.0),
        PVector(0.0, 0.0, 1.0)
    )

    # 面上非顶点 -> 在 circumsphere 内部
    @test insphere(t, PVector(0.5, 0.5, 0.0)) == Interior()  # z=0 面
    @test insphere(t, PVector(0.5, 0.0, 0.0)) == Interior()  # z=0 面
    @test insphere(t, PVector(0.0, 0.5, 0.0)) == Interior()  # z=0 面
    @test insphere(t, PVector(0.0, 0.5, 0.5)) == Interior()  # x=0 面
    @test insphere(t, PVector(0.0, 0.0, 0.5)) == Interior()  # x=0 面
    @test insphere(t, PVector(0.5, 0.0, 0.5)) == Interior()  # y=0 面

    # 顶点 -> OnEdge
    @test insphere(t, PVector(0.0, 0.0, 0.0)) == OnEdge()
    @test insphere(t, PVector(1.0, 0.0, 0.0)) == OnEdge()
    @test insphere(t, PVector(0.0, 1.0, 0.0)) == OnEdge()
    @test insphere(t, PVector(0.0, 0.0, 1.0)) == OnEdge()

    # 内部点
    @test insphere(t, PVector(0.1, 0.1, 0.1)) == Interior()
    @test insphere(t, PVector(0.01, 0.01, 0.01)) == Interior()

    # 外部点; (1,1,1) 距中心 sqrt(0.75) = 半径, 恰在 circumsphere 上 -> OnEdge
    @test insphere(t, PVector(1.0, 1.0, 1.0)) == OnEdge()
    @test insphere(t, PVector(2.0, 0.0, 0.0)) == Exterior()
    @test insphere(t, PVector(-0.1, -0.1, -0.1)) == Exterior()
end

@testset "Tetrahedron - Orientation" begin
    # 正方向
    t_pos = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(1.0, 0.0, 0.0),
        PVector(0.0, 1.0, 0.0),
        PVector(0.0, 0.0, 1.0)
    )
    @test orientation(t_pos) == PositivelyOriented()
    @test PhysicalMeshes.orient(t_pos) > 0

    # 负方向（交换两个顶点）
    t_neg = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(0.0, 1.0, 0.0),
        PVector(1.0, 0.0, 0.0),
        PVector(0.0, 0.0, 1.0)
    )
    @test orientation(t_neg) == NegativelyOriented()
    @test PhysicalMeshes.orient(t_neg) < 0
    
    # 负方向的体积应该与正方向相同（绝对值）
    @test PhysicalMeshes.volume(t_pos) == PhysicalMeshes.volume(t_neg)
end

@testset "Tetrahedron - Standard" begin
    # 单位立方体的一个角
    t_unit = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(1.0, 0.0, 0.0),
        PVector(0.0, 1.0, 0.0),
        PVector(0.0, 0.0, 1.0)
    )
    @test PhysicalMeshes.volume(t_unit) == 1.0 / 6.0
    @test centroid(t_unit) == PVector(0.25, 0.25, 0.25)

    # 验证 circumcenter
    center = circumcenter(t_unit)
    # 验证 center 到四个顶点距离相等
    distances = [norm(center - v) for v in [t_unit.a, t_unit.b, t_unit.c, t_unit.d]]
    @test all(d -> d ≈ distances[1], distances)
    
    # 验证 circumcenter 位置（对于单位直角四面体，circumcenter 应该在 (0.5, 0.5, 0.5)）
    @test center ≈ PVector(0.5, 0.5, 0.5)
end

@testset "Tetrahedron - With units" begin
    using PhysicalParticles.Unitful
    
    t_unitful = Tetrahedron(
        PVector(0.0, 0.0, 0.0, u"m"),
        PVector(1.0, 0.0, 0.0, u"m"),
        PVector(0.0, 1.0, 0.0, u"m"),
        PVector(0.0, 0.0, 1.0, u"m")
    )
    @test PhysicalMeshes.volume(t_unitful) == (1.0 / 6.0)u"m^3"
    @test centroid(t_unitful) == PVector(0.25, 0.25, 0.25, u"m")
end

@testset "Tetrahedron - Other edge cases" begin
    # 非常小的四面体
    t_small = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(1e-10, 0.0, 0.0),
        PVector(0.0, 1e-10, 0.0),
        PVector(0.0, 0.0, 1e-10)
    )
    @test PhysicalMeshes.volume(t_small) ≈ (1e-30) / 6.0
    @test orientation(t_small) == PositivelyOriented()
    
    # 非常大的四面体
    t_large = Tetrahedron(
        PVector(0.0, 0.0, 0.0),
        PVector(1e10, 0.0, 0.0),
        PVector(0.0, 1e10, 0.0),
        PVector(0.0, 0.0, 1e10)
    )
    @test PhysicalMeshes.volume(t_large) ≈ (1e30) / 6.0
    @test orientation(t_large) == PositivelyOriented()
end
