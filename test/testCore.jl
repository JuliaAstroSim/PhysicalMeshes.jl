using Test
using PhysicalParticles
using PhysicalMeshes: orient_kernel

@testset "core/orient.jl" begin
    @testset "2D orient_kernel: CCW vs CW vs collinear" begin
        # CCW triangle (counter-clockwise) → positive
        @test orient_kernel(
            PVector2D(0.0, 0.0),
            PVector2D(1.0, 0.0),
            PVector2D(0.0, 1.0),
        ) > 0
        # CW triangle → negative (sign flips vs CCW)
        @test orient_kernel(
            PVector2D(0.0, 0.0),
            PVector2D(0.0, 1.0),
            PVector2D(1.0, 0.0),
        ) < 0
        # Collinear → zero (kernel of signed area)
        @test isapprox(orient_kernel(
            PVector2D(0.0, 0.0),
            PVector2D(1.0, 1.0),
            PVector2D(2.0, 2.0),
        ), 0.0; atol=1e-12)
    end

    @testset "2D orient returns (0,0,z)" begin
        # The 2D `orient` returns a PVector with .z holding twice the signed area.
        z = orient(
            PVector2D(0.0, 0.0),
            PVector2D(2.0, 0.0),
            PVector2D(0.0, 3.0),
        )
        # Area = 0.5*2*3 = 3, so twice the area = 6
        @test z.x == 0
        @test z.y == 0
        @test z.z ≈ 6.0
    end

    @testset "3D orient: triangle normal direction" begin
        # CCW triangle in xy-plane (z=0), normal points to +z (right-hand rule)
        n_ccw = orient(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
        )
        @test n_ccw.z > 0
        @test isapprox(n_ccw.x, 0.0; atol=1e-12)
        @test isapprox(n_ccw.y, 0.0; atol=1e-12)

        # CW triangle in xy-plane → normal points to -z
        n_cw = orient(
            PVector(0.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
            PVector(1.0, 0.0, 0.0),
        )
        @test n_cw.z < 0
    end

    @testset "3D orient: collinear points → zero magnitude" begin
        n = orient(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(2.0, 0.0, 0.0),
        )
        @test isapprox(n.x, 0.0; atol=1e-12)
        @test isapprox(n.y, 0.0; atol=1e-12)
        @test isapprox(n.z, 0.0; atol=1e-12)
    end

    @testset "3D orient: magnitude equals 2x triangle area" begin
        # Right triangle with legs 3 (x) and 4 (y) in the xy-plane: area = 6
        n = orient(
            PVector(0.0, 0.0, 0.0),
            PVector(3.0, 0.0, 0.0),
            PVector(0.0, 4.0, 0.0),
        )
        @test isapprox(sqrt(n.x^2 + n.y^2 + n.z^2), 2.0 * 6.0; atol=1e-12)
    end

    @testset "4-point orient_kernel: signed tetrahedron volume" begin
        # Regular tetrahedron-ish: determinant ≠ 0
        v = orient_kernel(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
            PVector(0.0, 0.0, 1.0),
        )
        @test v != 0
        # Coplanar points → zero
        v_coplanar = orient_kernel(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
            PVector(0.5, 0.5, 0.0),
        )
        @test isapprox(v_coplanar, 0.0; atol=1e-12)
    end

    @testset "4-point orient: sign flips under point reversal" begin
        v1 = orient(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
            PVector(0.0, 0.0, 1.0),
        )
        # Swap two interior points → sign flips
        v2 = orient(
            PVector(0.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 0.0, 1.0),
        )
        @test v1 * v2 < 0
    end
end

@testset "core/circumcenter.jl" begin
    @testset "2D circumcenter: equilateral triangle" begin
        # Equilateral triangle with vertices on the unit circle:
        #   a = ( 1,   0), b = (-1/2,  √3/2), c = (-1/2, -√3/2)
        # The circumcenter is at the origin (0, 0) and the circumradius
        # (distance from cc to each vertex) is 1.
        a = PVector2D(1.0, 0.0)
        b = PVector2D(-0.5,  sqrt(3)/2)
        c = PVector2D(-0.5, -sqrt(3)/2)
        cc = circumcenter(a, b, c)
        @test isapprox(cc.x, 0.0; atol=1e-12)
        @test isapprox(cc.y, 0.0; atol=1e-12)
        # Distance from cc to each vertex equals the circumradius (= 1)
        @test isapprox(sqrt((cc.x - a.x)^2 + (cc.y - a.y)^2), 1.0; atol=1e-12)
        @test isapprox(sqrt((cc.x - b.x)^2 + (cc.y - b.y)^2), 1.0; atol=1e-12)
        @test isapprox(sqrt((cc.x - c.x)^2 + (cc.y - c.y)^2), 1.0; atol=1e-12)
    end

    @testset "2D circumcenter: right triangle" begin
        # Right triangle with legs on axes → circumcenter is at (0.5, 0.5)
        a = PVector2D(0.0, 0.0)
        b = PVector2D(1.0, 0.0)
        c = PVector2D(0.0, 1.0)
        cc = circumcenter(a, b, c)
        @test isapprox(cc.x, 0.5; atol=1e-12)
        @test isapprox(cc.y, 0.5; atol=1e-12)
    end

    @testset "2D circumcenter: collinear → ArgumentError" begin
        @test_throws ArgumentError circumcenter(
            PVector2D(0.0, 0.0),
            PVector2D(1.0, 0.0),
            PVector2D(2.0, 0.0),
        )
    end

    @testset "3D circumcenter: triangle in 3D" begin
        # Same triangle as the equilateral test, lifted into z=0 plane
        a = PVector(1.0, 0.0, 0.0)
        b = PVector(-0.5,  sqrt(3)/2, 0.0)
        c = PVector(-0.5, -sqrt(3)/2, 0.0)
        cc = circumcenter(a, b, c)
        @test isapprox(cc.x, 0.0; atol=1e-12)
        @test isapprox(cc.y, 0.0; atol=1e-12)
        @test isapprox(cc.z, 0.0; atol=1e-12)
    end

    @testset "3D circumcenter: collinear → ArgumentError" begin
        @test_throws ArgumentError circumcenter(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(2.0, 0.0, 0.0),
        )
    end

    @testset "3D circumcenter of tetrahedron (4 points)" begin
        # Regular tetrahedron with vertices at (±1, ±1, ±1) with even parity.
        # The circumcenter is at the origin and the circumradius is sqrt(3).
        a = PVector( 1.0,  1.0,  1.0)
        b = PVector( 1.0, -1.0, -1.0)
        c = PVector(-1.0,  1.0, -1.0)
        d = PVector(-1.0, -1.0,  1.0)
        cc = circumcenter(a, b, c, d)
        @test isapprox(cc.x, 0.0; atol=1e-12)
        @test isapprox(cc.y, 0.0; atol=1e-12)
        @test isapprox(cc.z, 0.0; atol=1e-12)
        # Distance from cc to any vertex = sqrt(3)
        @test isapprox(sqrt((cc.x - a.x)^2 + (cc.y - a.y)^2 + (cc.z - a.z)^2),
                       sqrt(3.0); atol=1e-12)
        @test isapprox(sqrt((cc.x - b.x)^2 + (cc.y - b.y)^2 + (cc.z - b.z)^2),
                       sqrt(3.0); atol=1e-12)
    end

    @testset "3D circumcenter of tetrahedron: coplanar → ArgumentError" begin
        # All four points in z=0 plane
        @test_throws ArgumentError circumcenter(
            PVector(0.0, 0.0, 0.0),
            PVector(1.0, 0.0, 0.0),
            PVector(0.0, 1.0, 0.0),
            PVector(0.5, 0.5, 0.0),
        )
    end
end

@testset "core/predicates.jl" begin
    @testset "incircle: d inside CCW triangle's circumcircle" begin
        # CCW triangle: a=(0,0), b=(1,0), c=(0,1) → circumcircle center (0.5, 0.5), R = sqrt(0.5)
        # Pick a point well inside the triangle: clearly inside the circumcircle.
        d_in = PVector2D(0.25, 0.25)
        @test incircle(
            PVector2D(0.0, 0.0),
            PVector2D(1.0, 0.0),
            PVector2D(0.0, 1.0),
            d_in,
        ) isa Interior
    end

    @testset "incircle: d far outside the circumcircle" begin
        # Same CCW triangle; pick d far away from origin → outside the circumcircle.
        d_out = PVector2D(10.0, 10.0)
        @test incircle(
            PVector2D(0.0, 0.0),
            PVector2D(1.0, 0.0),
            PVector2D(0.0, 1.0),
            d_out,
        ) isa Exterior
    end

    @testset "incircle: orientation-independent" begin
        # For CW triangle (a, c, b), Shewchuk's predicate should still
        # classify based on *absolute* position relative to the circumcircle.
        d_in = PVector2D(0.25, 0.25)
        d_out = PVector2D(10.0, 10.0)
        # Same triangle with reversed order (CW)
        @test incircle(
            PVector2D(0.0, 0.0),
            PVector2D(0.0, 1.0),
            PVector2D(1.0, 0.0),
            d_in,
        ) isa Interior
        @test incircle(
            PVector2D(0.0, 0.0),
            PVector2D(0.0, 1.0),
            PVector2D(1.0, 0.0),
            d_out,
        ) isa Exterior
    end

    @testset "insphere: returns a valid predicate result" begin
        # Regular tetrahedron with vertices at (±1, ±1, ±1) with even parity.
        a = PVector( 1.0,  1.0,  1.0)
        b = PVector( 1.0, -1.0, -1.0)
        c = PVector(-1.0,  1.0, -1.0)
        d = PVector(-1.0, -1.0,  1.0)
        # The predicate must return one of the three valid trait instances.
        # (Shewchuk's sign convention is orientation-dependent, so we don't
        # assert Interior/Exterior directly — instead we check that the
        # result is consistent and changes for far-away points.)
        e_center = PVector(0.0, 0.0, 0.0)
        r_center = insphere(a, b, c, d, e_center)
        @test r_center isa Union{Interior, Exterior, OnEdge}

        # Far-away point must give a different classification than the center
        # (origin is the circumsphere center → inside; 10,10,10 is far outside).
        e_far = PVector(10.0, 10.0, 10.0)
        r_far = insphere(a, b, c, d, e_far)
        @test r_far isa Union{Interior, Exterior, OnEdge}
        @test r_center != r_far
    end
end
