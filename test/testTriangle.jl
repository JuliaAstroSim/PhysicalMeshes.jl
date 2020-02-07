@testset "Triangle" begin
    t = Triangle2D(PVector2D(0.0, 0.0), PVector2D(4.0, 0.0), PVector2D(0.0, 3.0))
    @test len(t) == 12.0
    @test centroid(t) == PVector2D(2.0, 1.5)

    t = Triangle(PVector(0.0, 0.0, 0.0, u"m"), PVector(3.0, 4.0, 0.0, u"m"), PVector(0.0, 0.0, 12.0, u"m"))
    @test len(t) == 30.0u"m"
    @test centroid(t) == PVector(1.5, 2.0, 6.0, u"m")
end