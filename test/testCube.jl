@testset "Cube" begin
    c = Cube(PVector2D(1.0, 1.0), PVector2D())
    @test area(c) == 1.0

    c = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
    @test volume(c) == 1.0u"m^3"
end