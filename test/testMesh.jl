@testset "Mesh" begin
    # data
    pos = [
        PVector(-1.0, -1.0, -1.0),
        PVector(-1.0, +1.0, -1.0),
        PVector(+1.0, -1.0, -1.0),
        PVector(+1.0, +1.0, -1.0),
        PVector(-1.0, -1.0, +1.0),
        PVector(-1.0, +1.0, +1.0),
        PVector(+1.0, -1.0, +1.0),
        PVector(+1.0, +1.0, +1.0),
    ]
    dataArray = [Ball() for i in 1:8]
    assign_particles(dataArray, :Pos, pos)
    assign_particles(dataArray, :Mass, 1.0)

    @testset "Static Cartesian Mesh" begin
        m = MeshCartesianStatic(dataArray)

        # Check mass assignment
        @test sum(m.rho) * m.config.Δ[1]^3 ≈ 8.0
    end
end