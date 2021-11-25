@testset "Cube" begin
    c = Cube(PVector2D(1.0, 1.0), PVector2D())
    @test PhysicalMeshes.area(c) == 1.0

    c = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
    @test PhysicalMeshes.volume(c) == 1.0u"m^3"
end