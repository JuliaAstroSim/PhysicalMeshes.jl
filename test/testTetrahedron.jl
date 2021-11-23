@testset "Tetrahedron" begin
    t = Tetrahedron(
        PVector(1.0, 0.0, 1.0),
        PVector(1.0, 1.0, 0.0),
        PVector(0.0, 1.0, 1.0),
        PVector(1.0, 1.0, 1.0)
    )

    @test centroid(t) == PVector(3.0, 3.0, 3.0) / 4.0
    @test circumcenter(t) == PVector(0.5, 0.5, 0.5)
    
    @test insphere(t, PVector(0.5, 0.5, 0.5)) == Interior()
    @test insphere(t, PVector()) == OnEdge()
    @test insphere(t, PVector(2.0, 0.0, 0.0)) == Exterior()

    @test volume(t) == 1.0 / 6.0

    @test orientation(t) == PositivelyOriented()
end