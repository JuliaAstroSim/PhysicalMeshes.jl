using Test
using PhysicalParticles

@testset "Plane" begin
    # Standard XY plane (z=0), defined by three points: origin, x̂, ŷ
    plane_xy = Plane(
        PVector(0.0, 0.0, 0.0),
        PVector(1.0, 0.0, 0.0),
        PVector(0.0, 1.0, 0.0),
    )

    @testset "normal of canonical XY plane" begin
        n = PhysicalMeshes.normal(plane_xy)
        # The XY plane normal is ±ẑ; magnitude must be ~1 and z component dominant
        @test isapprox(hypot(n.x, n.y, n.z), 1.0; atol=1e-12)
        @test isapprox(abs(n.z), 1.0; atol=1e-12)
        @test isapprox(n.x, 0.0; atol=1e-12)
        @test isapprox(n.y, 0.0; atol=1e-12)
    end

    @testset "normal of XZ plane (y=0)" begin
        plane_xz = Plane(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 0.0, 1.0),
        )
        n = PhysicalMeshes.normal(plane_xz)
        @test isapprox(hypot(n.x, n.y, n.z), 1.0; atol=1e-12)
        @test isapprox(abs(n.y), 1.0; atol=1e-12)
    end

    @testset "normal of tilted plane (non-canonical)" begin
        # Plane through (0,0,0), (1,1,0), (0,1,1) → normal ∝ (-1, 1, 1)
        plane = Plane(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 1.0, 0.0),
            PVector(0.0, 1.0, 1.0),
        )
        n = PhysicalMeshes.normal(plane)
        @test isapprox(hypot(n.x, n.y, n.z), 1.0; atol=1e-12)
        # Each component should have equal magnitude
        @test isapprox(abs(n.x), abs(n.y); atol=1e-12)
        @test isapprox(abs(n.y), abs(n.z); atol=1e-12)
    end

    @testset "normal is right-hand-rule oriented" begin
        # Right-hand rule: cross(b-a, c-a) for plane (a, b, c).
        # Take plane (0,0,0), (1,0,0), (0,1,0) → cross = (0,0,1) → +ẑ.
        plane = Plane(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
        )
        n = PhysicalMeshes.normal(plane)
        @test n.z > 0  # +ẑ orientation
    end

    @testset "distance to plane" begin
        plane = Plane(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
        )
        # `distance` is imported from PhysicalParticles in PhysicalMeshes.jl,
        # so the method `distance(::AbstractPoint3D, ::Plane)` extends the
        # global `distance` function. Unqualified call dispatches correctly
        # after `using PhysicalMeshes` + `using PhysicalParticles`.
        @test isapprox(distance(PVector(0.5, 0.5, 0.0), plane), 0.0; atol=1e-12)
        @test isapprox(distance(PVector(0.0, 0.0, 0.0), plane), 0.0; atol=1e-12)
        @test isapprox(distance(PVector(0.0, 0.0, 1.0), plane), 1.0; atol=1e-12)
        @test isapprox(distance(PVector(1.0, 1.0, -2.5), plane), 2.5; atol=1e-12)
        # Symmetric under sign flip of the perpendicular component
        @test isapprox(
            distance(PVector(0.0, 0.0,  3.0), plane),
            distance(PVector(0.0, 0.0, -3.0), plane);
            atol=1e-12,
        )
    end

    @testset "coplanar predicate" begin
        plane = Plane(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
        )
        # Strict coplanar points
        @test coplanar(PVector(0.5, 0.5, 0.0), plane, 1e-12)
        # Within threshold
        @test coplanar(PVector(0.5, 0.5, 1e-8), plane, 1e-6)
        # Outside threshold
        @test !coplanar(PVector(0.5, 0.5, 1.0), plane, 0.1)
        @test !coplanar(PVector(0.5, 0.5, 1.0), plane, 0.5)
    end

    @testset "Plane struct fields" begin
        a = PVector(1.0, 2.0, 3.0)
        b = PVector(4.0, 5.0, 6.0)
        c = PVector(7.0, 8.0, 10.0)
        plane = Plane(a, b, c)
        @test plane.a == a
        @test plane.b == b
        @test plane.c == c
    end
end
