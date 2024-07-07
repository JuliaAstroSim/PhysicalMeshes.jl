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
    @test_broken is_inbound(PVector(), polygon_nonconvex)
    @test_broken !is_inbound(PVector(2.0, 0.0, 0.0), polygon_nonconvex)
end

@testset "Polygon3D Unitful" begin
    polygon = Polygon3D([PVector(1.0, 1.0, 0.0, u"m"), PVector(1.0, -1.0, 0.0, u"m"), PVector(-1.0, -1.0, 0.0, u"m"), PVector(-1.0, 1.0, 0.0, u"m")])
    @test isconvex(polygon)
    @test is_inbound(PVector(u"m"), polygon)
    @test !is_inbound(PVector(2.0, 0.0, 0.0, u"m"), polygon)

    polygon_nonconvex = Polygon3D([PVector(1.0, 1.0, 0.0, u"m"), PVector(0.5, 0.0, 0.0, u"m"), PVector(1.0, -1.0, 0.0, u"m"), PVector(-1.0, -1.0, 0.0, u"m"), PVector(-1.0, 1.0, 0.0, u"m")])
    @test !isconvex(polygon_nonconvex)
    @test_broken is_inbound(PVector(u"m"), polygon_nonconvex)
    @test_broken !is_inbound(PVector(2.0, 0.0, 0.0, u"m"), polygon_nonconvex)
end