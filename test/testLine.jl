@testset "Line" begin
    @testset "Line2D" begin
        line = Line2D(PVector2D(), PVector2D(1.0, 0.0))
        @test len(line) == 1.0
        @test midpoint(line) == PVector2D(0.5, 0.0)

        m = Line2D(PVector2D(2.0, 2.0, u"m"), PVector2D(2.0, 2.0, u"m"))
        n = Line2D(PVector2D(1.0, 1.0, u"m"), PVector2D(1.0, 1.0, u"m"))

        @test m + n == Line(PVector2D(3.0, 3.0, u"m"), PVector2D(3.0, 3.0, u"m"))
        @test m - n == Line(PVector2D(1.0, 1.0, u"m"), PVector2D(1.0, 1.0, u"m"))
        @test n * 2 == m
        @test 2 * n == m

        @test n + PVector2D(1.0, 1.0, u"m") == m
        @test m - PVector2D(1.0, 1.0, u"m") == n
        @test m / 2 == n
    end
    
    @testset "Line" begin
        line = Line(PVector(u"m"), PVector(1.0, 0.0, 0.0, u"m"))
        @test len(line) == 1.0u"m"
        @test midpoint(line) == PVector(0.5, 0.0, 0.0, u"m")

        m = Line(PVector(2.0, 2.0, 2.0), PVector(2.0, 2.0, 2.0))
        n = Line(PVector(1.0, 1.0, 1.0), PVector(1.0, 1.0, 1.0))

        @test m + n == Line(PVector(3.0, 3.0, 3.0), PVector(3.0, 3.0, 3.0))
        @test m - n == Line(PVector(1.0, 1.0, 1.0), PVector(1.0, 1.0, 1.0))
        @test n * 2 == m
        @test 2 * n == m

        @test n + PVector(1.0, 1.0, 1.0) == m
        @test m - PVector(1.0, 1.0, 1.0) == n
        @test m / 2 == n
    end
end