@testset "Line" begin
    line = Line2D(PVector2D(), PVector2D(1.0, 0.0))
    @test len(line) == 1.0

    line = Line(PVector(u"m"), PVector(1.0, 0.0, 0.0, u"m"))
    @test len(line) == 1.0u"m"
end