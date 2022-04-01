@testset "Cube" begin
    c = Cube(PVector2D(1.0, 1.0), PVector2D())
    @test PhysicalMeshes.area(c) == 1.0

    c = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
    @test PhysicalMeshes.volume(c) == 1.0u"m^3"

    @test interior(c, PVector(0.5, 0.5, 0.5, u"m"))
    @test exterior(c, PVector(1.5, 1.5, 1.5, u"m"))
    @test !interior(c, PVector(u"m"))
    @test !exterior(c, PVector(u"m"))
end