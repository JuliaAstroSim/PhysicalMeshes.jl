@testset "Sphere" begin
    s = Sphere(PVector(), 1.0)
    @test interior(s, PVector(0.1, 0.1, 0.1))
    @test exterior(s, PVector(1.0, 1.0, 1.0))
    @test !interior(s, PVector(1.0, 0.0, 0.0))
    @test !exterior(s, PVector(1.0, 0.0, 0.0))
end