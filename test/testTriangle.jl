@testset "Triangle" begin
    @testset "Triangle2D" begin    
        t = Triangle2D(PVector2D(4.0, 3.0), PVector2D(4.0, 0.0), PVector2D(0.0, 3.0))

        @test orientation(t) == NegativelyOriented()

        @test len(t) == 12.0
        @test area(t) == 6.0

        @test centroid(t) == PVector2D(8.0/3.0, 2.0)
        @test circumcenter(t) == PVector2D(2.0, 1.5)

        @test incircle(t, PVector2D(3.0, 2.0)) == Inner()
        @test incircle(t, PVector2D()) == OnEdge()
        @test incircle(t, PVector2D(6.0, 0.0)) == Outter()
    end

    @testset "Triangle" begin
        t = Triangle(PVector(0.0, 0.0, 0.0, u"m"), PVector(3.0, 4.0, 0.0, u"m"), PVector(3.0, 4.0, 12.0, u"m"))

        @test len(t) == 30.0u"m"
        @test area(t) == 30.0u"m^2"

        @test centroid(t) == PVector(2.0, 8.0 / 3.0, 4.0, u"m")
        @test circumcenter(t) == PVector(1.5, 2.0, 6.0, u"m")
    end
end