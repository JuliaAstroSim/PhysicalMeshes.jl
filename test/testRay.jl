@testset "Ray2D Unitless" begin
    ray = Ray2D(PVector2D(), PVector2D(1.0, 1.0))
    @test ray.theta == π/4

    line = Line2D(PVector2D(2.0, 1.0), PVector2D(2.0, 3.0))
    hit, intersection = PhysicalMeshes.intersect(ray, line)
    @test intersection == PVector(2.0, 2.0)

    N = normal(line)
    @test N == PVector(-1.0, 0.0)

    ray_reflect = reflect(ray, line)
    @test ray_reflect.x == PVector(2.0, 2.0)
    @test ray_reflect.theta == 3/4 * π

    @test isnothing(reflect(ray, Line2D(PVector(2, 1), PVector(2, 0))))
end

@testset "Ray2D Unitful" begin
    ray = Ray2D(PVector2D(u"m"), PVector2D(1.0u"m", 1.0u"m"))
    @test ray.theta == π/4

    line = Line2D(PVector2D(2.0u"m", 1.0u"m"), PVector2D(2.0u"m", 3.0u"m"))
    hit, intersection = PhysicalMeshes.intersect(ray, line)
    @test intersection == PVector(2.0u"m", 2.0u"m")

    N = normal(line)
    @test N == PVector(-1.0u"m", 0.0u"m")

    ray_reflect = reflect(ray, line)
    @test ray_reflect.x == PVector(2.0u"m", 2.0u"m")
    @test ray_reflect.theta == 3/4 * π

    @test isnothing(reflect(ray, Line2D(PVector(2u"m", 1u"m"), PVector(2u"m", 0u"m"))))
end

@testset "Ray3D Unitless" begin
    ray = Ray3D(PVector(), PVector(0.0, 1.0, 1.0))
    plane = Plane(PVector(1, 1, 1), PVector(1, 1, -1), PVector(-1, 1, 1))
    @test normal(plane) == PVector(0.0, 1.0, 0.0)
    hit, intersection_point = intersect(ray, plane)
    @test intersection_point == PVector(0.0, 1.0, 1.0)
    ray_reflect = reflect(ray, plane)
    @test ray_reflect.n == PVector(0.0, -1.0, 1.0)
    
    #: zero incident angle
    ray = Ray3D(PVector(), PVector(0.0, 1.0, 0.0))
    plane = Plane(PVector(1, 1, 1), PVector(1, 1, -1), PVector(-1, 1, 1))
    @test normal(plane) == PVector(0.0, 1.0, 0.0)
    hit, intersection_point = intersect(ray, plane)
    @test intersection_point == PVector(0.0, 1.0, 0.0)
    ray_reflect = reflect(ray, plane)
    @test ray_reflect.n == PVector(0.0, -1.0, 0.0)
end

@testset "Ray3D Unitful" begin
    ray = Ray3D(PVector(u"m"), PVector(0.0, 1.0, 1.0, u"m"))
    plane = Plane(PVector(1, 1, 1, u"m"), PVector(1, 1, -1, u"m"), PVector(-1, 1, 1, u"m"))
    @test normal(plane) == PVector(0.0, 1.0, 0.0, u"m")
    hit, intersection_point = intersect(ray, plane)
    @test intersection_point == PVector(0.0, 1.0, 1.0, u"m")
    ray_reflect = reflect(ray, plane)
    @test ray_reflect.n == PVector(0.0, -1.0, 1.0, u"m")
    
    #: zero incident angle
    ray = Ray3D(PVector(u"m"), PVector(0.0, 1.0, 0.0, u"m"))
    plane = Plane(PVector(1, 1, 1, u"m"), PVector(1, 1, -1, u"m"), PVector(-1, 1, 1, u"m"))
    @test normal(plane) == PVector(0.0, 1.0, 0.0, u"m")
    hit, intersection_point = intersect(ray, plane)
    @test intersection_point == PVector(0.0, 1.0, 0.0, u"m")
    ray_reflect = reflect(ray, plane)
    @test ray_reflect.n == PVector(0.0, -1.0, 0.0, u"m")
end