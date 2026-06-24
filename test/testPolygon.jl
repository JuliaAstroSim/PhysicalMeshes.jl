@testset "Polygon2D Unitless" begin
    polygon = Polygon2D([PVector2D(1.0, 1.0), PVector2D(1.0, -1.0), PVector2D(-1.0, -1.0), PVector2D(-1.0, 1.0)])
    @test isconvex(polygon)
    @test is_inbound(PVector2D(), polygon)
    @test !is_inbound(PVector2D(2.0, 0.0), polygon)

    polygon_nonconvex = Polygon2D([PVector2D(1.0, 1.0), PVector2D(0.5, 0.0), PVector2D(1.0, -1.0), PVector2D(-1.0, -1.0), PVector2D(-1.0, 1.0)])
    @test !isconvex(polygon_nonconvex)
    @test is_inbound(PVector2D(), polygon_nonconvex)
    @test !is_inbound(PVector2D(2.0, 0.0), polygon_nonconvex)
end

@testset "Polygon2D Unitful" begin
    polygon = Polygon2D([PVector2D(1.0, 1.0, u"m"), PVector2D(1.0, -1.0, u"m"), PVector2D(-1.0, -1.0, u"m"), PVector2D(-1.0, 1.0, u"m")])
    @test isconvex(polygon)
    @test is_inbound(PVector2D(u"m"), polygon)
    @test !is_inbound(PVector2D(2.0, 0.0, u"m"), polygon)

    polygon_nonconvex = Polygon2D([PVector2D(1.0, 1.0, u"m"), PVector2D(0.5, 0.0, u"m"), PVector2D(1.0, -1.0, u"m"), PVector2D(-1.0, -1.0, u"m"), PVector2D(-1.0, 1.0, u"m")])
    @test !isconvex(polygon_nonconvex)
    @test is_inbound(PVector2D(u"m"), polygon_nonconvex)
    @test !is_inbound(PVector2D(2.0, 0.0, u"m"), polygon_nonconvex)
end

@testset "Polygon3D Unitless" begin
    polygon = Polygon3D([PVector(1.0, 1.0, 0.0), PVector(1.0, -1.0, 0.0), PVector(-1.0, -1.0, 0.0), PVector(-1.0, 1.0, 0.0)])
    @test isconvex(polygon)
    @test is_inbound(PVector(), polygon)
    @test !is_inbound(PVector(2.0, 0.0, 0.0), polygon)

    polygon_nonconvex = Polygon3D([PVector(1.0, 1.0, 0.0), PVector(0.5, 0.0, 0.0), PVector(1.0, -1.0, 0.0), PVector(-1.0, -1.0, 0.0), PVector(-1.0, 1.0, 0.0)])
    @test !isconvex(polygon_nonconvex)
    @test is_inbound(PVector(), polygon_nonconvex)
    @test !is_inbound(PVector(2.0, 0.0, 0.0), polygon_nonconvex)
end

@testset "Polygon3D Unitful" begin
    polygon = Polygon3D([PVector(1.0, 1.0, 0.0, u"m"), PVector(1.0, -1.0, 0.0, u"m"), PVector(-1.0, -1.0, 0.0, u"m"), PVector(-1.0, 1.0, 0.0, u"m")])
    @test isconvex(polygon)
    @test is_inbound(PVector(u"m"), polygon)
    @test !is_inbound(PVector(2.0, 0.0, 0.0, u"m"), polygon)

    polygon_nonconvex = Polygon3D([PVector(1.0, 1.0, 0.0, u"m"), PVector(0.5, 0.0, 0.0, u"m"), PVector(1.0, -1.0, 0.0, u"m"), PVector(-1.0, -1.0, 0.0, u"m"), PVector(-1.0, 1.0, 0.0, u"m")])
    @test !isconvex(polygon_nonconvex)
    @test is_inbound(PVector(u"m"), polygon_nonconvex)
    @test !is_inbound(PVector(2.0, 0.0, 0.0, u"m"), polygon_nonconvex)
end

# 凸多边形测试 - 2D
@testset "Polygon2D Convex" begin
    # 正方形
    square = Polygon2D([PVector2D(1.0, 1.0), PVector2D(1.0, -1.0), PVector2D(-1.0, -1.0), PVector2D(-1.0, 1.0)])
    @test isconvex(square)
    @test is_inbound(PVector2D(0.0, 0.0), square)
    @test is_inbound(PVector2D(0.5, 0.5), square)
    @test !is_inbound(PVector2D(2.0, 0.0), square)
    @test !is_inbound(PVector2D(0.0, 2.0), square)

    # 三角形
    triangle = Polygon2D([PVector2D(0.0, 0.0), PVector2D(2.0, 0.0), PVector2D(0.0, 2.0)])
    @test isconvex(triangle)
    @test is_inbound(PVector2D(0.5, 0.5), triangle)
    @test !is_inbound(PVector2D(2.0, 2.0), triangle)

    # 五边形
    pentagon = polygon_regular(5, 1.0)
    @test isconvex(pentagon)
    @test is_inbound(PVector2D(0.0, 0.0), pentagon)
end

# 非凸多边形测试 - 2D
@testset "Polygon2D Non-convex" begin
    # L 形
    l_shape = Polygon2D([
        PVector2D(0.0, 0.0), PVector2D(2.0, 0.0), PVector2D(2.0, 1.0),
        PVector2D(1.0, 1.0), PVector2D(1.0, 2.0), PVector2D(0.0, 2.0)
    ])
    @test !isconvex(l_shape)
    @test is_inbound(PVector2D(0.5, 0.5), l_shape)
    @test is_inbound(PVector2D(0.5, 1.5), l_shape)
    @test !is_inbound(PVector2D(1.5, 1.5), l_shape)  # 凹进去的部分

    # 星形（简化版）
    star = Polygon2D([
        PVector2D(0.0, 1.0), PVector2D(0.2, 0.2), PVector2D(1.0, 0.0),
        PVector2D(0.2, -0.2), PVector2D(0.0, -1.0), PVector2D(-0.2, -0.2),
        PVector2D(-1.0, 0.0), PVector2D(-0.2, 0.2)
    ])
    @test !isconvex(star)
    @test is_inbound(PVector2D(0.0, 0.0), star)
    @test !is_inbound(PVector2D(0.5, 0.5), star)  # 星形的凹处
end

# 凸多边形测试 - 3D
@testset "Polygon3D Convex" begin
    # xy 平面上的正方形
    square_3d = Polygon3D([
        PVector(1.0, 1.0, 0.0), PVector(1.0, -1.0, 0.0),
        PVector(-1.0, -1.0, 0.0), PVector(-1.0, 1.0, 0.0)
    ])
    @test isconvex(square_3d)
    @test is_inbound(PVector(0.0, 0.0, 0.0), square_3d)
    @test is_inbound(PVector(0.5, 0.5, 0.0), square_3d)
    @test !is_inbound(PVector(2.0, 0.0, 0.0), square_3d)
    @test !is_inbound(PVector(0.0, 0.0, 1.0), square_3d)  # 不在平面上

    # 斜平面上的多边形
    slanted = Polygon3D([
        PVector(1.0, 0.0, 1.0), PVector(0.0, 1.0, 1.0),
        PVector(-1.0, 0.0, 1.0), PVector(0.0, -1.0, 1.0)
    ])
    @test isconvex(slanted)
    @test is_inbound(PVector(0.0, 0.0, 1.0), slanted)
end

# 非凸多边形测试 - 3D（修复后的测试）
@testset "Polygon3D Non-convex" begin
    polygon_nonconvex = Polygon3D([
        PVector(1.0, 1.0, 0.0), PVector(0.5, 0.0, 0.0),
        PVector(1.0, -1.0, 0.0), PVector(-1.0, -1.0, 0.0),
        PVector(-1.0, 1.0, 0.0)
    ])
    @test !isconvex(polygon_nonconvex)

    # 这些测试当前标记为 @test_broken，需要修复后改为 @test
    @test is_inbound(PVector(0.0, 0.0, 0.0), polygon_nonconvex)
    @test !is_inbound(PVector(2.0, 0.0, 0.0), polygon_nonconvex)
    @test !is_inbound(PVector(0.0, 0.0, 1.0), polygon_nonconvex)  # 不在平面上

    # 测试凹进去的点
    @test !is_inbound(PVector(0.7, 0.0, 0.0), polygon_nonconvex)
end

# 边界测试
@testset "Polygon Boundary" begin
    square = Polygon2D([PVector2D(1.0, 1.0), PVector2D(1.0, -1.0), PVector2D(-1.0, -1.0), PVector2D(-1.0, 1.0)])

    # 恰好在边上
    @test is_inbound(PVector2D(1.0, 0.0), square)  # 或根据实现返回 false
    @test is_inbound(PVector2D(0.0, 1.0), square)

    # 恰好在顶点
    @test is_inbound(PVector2D(1.0, 1.0), square)
end

# polygon_regular 测试
@testset "polygon_regular" begin
    # 三角形
    tri = polygon_regular(3, 1.0)
    @test length(tri.points) == 3
    @test isconvex(tri)

    # 六边形
    hex = polygon_regular(6, 2.0)
    @test length(hex.points) == 6
    @test isconvex(hex)

    # 圆形近似（多边形）
    circle_approx = polygon_regular(100, 1.0)
    @test length(circle_approx.points) == 100
end

# polygon_rect 测试
@testset "polygon_rect" begin
    rect = polygon_rect(2.0, 3.0)
    @test length(rect.points) == 4
    @test isconvex(rect)
    @test is_inbound(PVector2D(0.0, 0.0), rect)
    @test is_inbound(PVector2D(0.9, 1.4), rect)
    @test !is_inbound(PVector2D(1.1, 0.0), rect)
end
