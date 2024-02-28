@testset "Mesh" begin
    @testset "Tools" begin
        @test zoom([1],2) == ones(Int, 2)
        @test zoom(ones(Int, 1, 1),2) == ones(Int, 2, 2)
        @test zoom(ones(Int, 1, 1, 1),2) == ones(Int, 2, 2, 2)

        @test zoom([1,2],2,interp=LinearInterpolation) ≈ [1.0, 4/3, 5/3, 2.0]

        a1 = [1,2,3,4]
        a2 = a1*a1'
        a3 = ones(Int,4,4,4) .* reshape(a1,1,1,4)

        # make sure that no error is thrown
        zoom([1 2; 2 3],2,interp=LinearInterpolation)
        zoom([1 2; 2 3],2,interp=LinearInterpolation)

        @test shrink(a1,2) == [2, 4]
        @test shrink(a2,2) == [4 8; 8 16]
        @test shrink(a3,2) == [2 2; 2 2;;; 4 4; 4 4;;;]
    end

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