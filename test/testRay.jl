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

# ========== 边界测试 ==========

@testset "Ray2D - Parallel reflection" begin
    # Ray 平行于反射面 - 水平线
    ray_parallel_h = Ray2D(PVector2D(0.0, 0.0), PVector2D(1.0, 0.0))  # 水平向右
    line_horizontal = Line2D(PVector2D(-1.0, 1.0), PVector2D(1.0, 1.0))  # 水平线 y=1
    hit, intersection = PhysicalMeshes.intersect(ray_parallel_h, line_horizontal)
    @test !hit  # 平行，无交点
    @test isnothing(reflect(ray_parallel_h, line_horizontal))  # 反射返回 nothing

    # Ray 平行于反射面 - 垂直线
    ray_parallel_v = Ray2D(PVector2D(0.0, 0.0), PVector2D(0.0, 1.0))  # 垂直向上
    line_vertical = Line2D(PVector2D(1.0, -1.0), PVector2D(1.0, 1.0))  # 垂直线 x=1
    hit, intersection = PhysicalMeshes.intersect(ray_parallel_v, line_vertical)
    @test !hit  # 平行，无交点
    @test isnothing(reflect(ray_parallel_v, line_vertical))  # 反射返回 nothing
end

@testset "Ray2D - Miss cases" begin
    # Ray 与线段不相交（在左侧）
    ray_right = Ray2D(PVector2D(0.0, 0.0), PVector2D(1.0, 0.0))  # 向右
    line_left = Line2D(PVector2D(-2.0, 1.0), PVector2D(-1.0, 2.0))  # 在左侧
    hit, intersection = PhysicalMeshes.intersect(ray_right, line_left)
    @test !hit
    @test isnothing(reflect(ray_right, line_left))

    # Ray 方向相反（向后）
    ray_left = Ray2D(PVector2D(0.0, 0.0), PVector2D(-1.0, 0.0))  # 向左
    line_right = Line2D(PVector2D(2.0, -1.0), PVector2D(2.0, 1.0))  # 在右侧
    hit, intersection = PhysicalMeshes.intersect(ray_left, line_right)
    @test !hit  # 向左的射线不会击中右侧的线
    @test isnothing(reflect(ray_left, line_right))

    # Ray 与线段延长线相交但不在线段范围内
    ray_up = Ray2D(PVector2D(0.0, 0.0), PVector2D(0.0, 1.0))  # 向上
    line_short = Line2D(PVector2D(1.0, 2.0), PVector2D(2.0, 2.0))  # 在右侧的水平短线
    hit, intersection = PhysicalMeshes.intersect(ray_up, line_short)
    @test !hit  # 向上射线不会击中右侧的线段
    @test isnothing(reflect(ray_up, line_short))
end

@testset "Ray2D - Zero direction" begin
    # 零方向向量
    ray_zero = Ray2D(PVector2D(0.0, 0.0), PVector2D(0.0, 0.0))
    @test iszero(ray_zero)
    @test ray_zero.theta == 0.0  # atan(0, 0) = 0

    line = Line2D(PVector2D(1.0, 1.0), PVector2D(2.0, 2.0))
    hit, intersection = PhysicalMeshes.intersect(ray_zero, line)
    # 零方向向量时，dot_ray_norm_AB 为 0，应该返回 false
    @test !hit
    @test isnothing(reflect(ray_zero, line))
end

@testset "Ray2D - Various reflection angles" begin
    # 45度入射到垂直线
    ray_45 = Ray2D(PVector2D(0.0, 0.0), PVector2D(1.0, 1.0))
    line_vertical = Line2D(PVector2D(2.0, 1.0), PVector2D(2.0, 3.0))  # 垂直线 x=2, 包含 y=2
    hit, intersection = PhysicalMeshes.intersect(ray_45, line_vertical)
    @test hit
    @test intersection == PVector2D(2.0, 2.0)
    ray_reflected = reflect(ray_45, line_vertical)
    @test ray_reflected.theta ≈ 3π/4  # 反射后向左上，135度

    # 垂直入射到水平线
    ray_vertical = Ray2D(PVector2D(0.0, 0.0), PVector2D(0.0, 1.0))  # 向上
    line_h = Line2D(PVector2D(-1.0, 2.0), PVector2D(1.0, 2.0))  # 水平线 y=2
    hit, intersection = PhysicalMeshes.intersect(ray_vertical, line_h)
    @test hit
    @test intersection == PVector2D(0.0, 2.0)
    ray_reflected_v = reflect(ray_vertical, line_h)
    @test ray_reflected_v.theta ≈ -π/2  # 反射后向下，-90度

    # 水平入射到垂直线
    ray_horizontal = Ray2D(PVector2D(0.0, 0.0), PVector2D(1.0, 0.0))  # 向右
    line_v = Line2D(PVector2D(2.0, -1.0), PVector2D(2.0, 1.0))  # 垂直线 x=2
    hit, intersection = PhysicalMeshes.intersect(ray_horizontal, line_v)
    @test hit
    @test intersection == PVector2D(2.0, 0.0)
    ray_reflected_h = reflect(ray_horizontal, line_v)
    @test ray_reflected_h.theta ≈ π  # 反射后向左，180度
end

@testset "Ray3D - Parallel reflection" begin
    # Ray 平行于平面 - x方向射线平行于xy平面
    ray_parallel_xy = Ray3D(PVector(0.0, 0.0, 0.0), PVector(1.0, 0.0, 0.0))  # x方向
    plane_parallel = Plane(PVector(0.0, 0.0, 1.0), PVector(1.0, 0.0, 1.0), PVector(0.0, 1.0, 1.0))  # z=1平面
    hit, point = PhysicalMeshes.intersect(ray_parallel_xy, plane_parallel)
    @test !hit  # 平行于平面，无交点
    @test isnothing(reflect(ray_parallel_xy, plane_parallel))

    # Ray 在平面内
    ray_in_plane = Ray3D(PVector(0.0, 0.0, 1.0), PVector(1.0, 0.0, 0.0))  # 起点在z=1平面内，x方向
    hit, point = PhysicalMeshes.intersect(ray_in_plane, plane_parallel)
    @test !hit  # 在平面内也视为平行，无交点
    @test isnothing(reflect(ray_in_plane, plane_parallel))
end

@testset "Ray3D - Miss cases" begin
    # Ray 与平面不相交（反向）
    ray_away = Ray3D(PVector(0.0, 0.0, 0.0), PVector(0.0, 0.0, -1.0))  # 向下
    plane_above = Plane(PVector(1.0, 1.0, 1.0), PVector(1.0, -1.0, 1.0), PVector(-1.0, 1.0, 1.0))  # z=1平面
    hit, point = PhysicalMeshes.intersect(ray_away, plane_above)
    @test !hit  # 向下射线不会击中上方的平面
    @test isnothing(reflect(ray_away, plane_above))

    # Ray 起点在平面后方
    ray_from_behind = Ray3D(PVector(0.0, 0.0, 2.0), PVector(0.0, 0.0, 1.0))  # 从z=2向上
    hit, point = PhysicalMeshes.intersect(ray_from_behind, plane_above)
    @test !hit  # 远离平面
    @test isnothing(reflect(ray_from_behind, plane_above))
end

@testset "Ray3D - Zero direction" begin
    # 零方向向量
    ray_zero = Ray3D(PVector(0.0, 0.0, 0.0), PVector(0.0, 0.0, 0.0))
    @test iszero(ray_zero)

    plane = Plane(PVector(1.0, 0.0, 0.0), PVector(0.0, 1.0, 0.0), PVector(0.0, 0.0, 1.0))
    hit, point = PhysicalMeshes.intersect(ray_zero, plane)
    # 零方向向量时，denom = 0，应该返回 false
    @test !hit
    @test isnothing(reflect(ray_zero, plane))
end

@testset "Ray3D - Various reflection angles" begin
    # 45度入射到xy平面
    ray_45 = Ray3D(PVector(0.0, 0.0, 0.0), PVector(0.0, 1.0, 1.0))
    plane_xy = Plane(PVector(1.0, 1.0, 1.0), PVector(1.0, -1.0, 1.0), PVector(-1.0, 1.0, 1.0))  # z=1平面
    hit, point = PhysicalMeshes.intersect(ray_45, plane_xy)
    @test hit
    @test point == PVector(0.0, 1.0, 1.0)
    ray_reflected = reflect(ray_45, plane_xy)
    @test ray_reflected.n == PVector(0.0, 1.0, -1.0)  # z方向反转（法向量方向）

    # 垂直入射（法向量方向）
    ray_normal = Ray3D(PVector(0.0, 0.0, 0.0), PVector(0.0, 0.0, 1.0))
    plane_z = Plane(PVector(1.0, 0.0, 1.0), PVector(0.0, 1.0, 1.0), PVector(-1.0, 0.0, 1.0))  # z=1平面
    hit, point = PhysicalMeshes.intersect(ray_normal, plane_z)
    @test hit
    @test point == PVector(0.0, 0.0, 1.0)
    ray_reflected_normal = reflect(ray_normal, plane_z)
    @test ray_reflected_normal.n == PVector(0.0, 0.0, -1.0)  # 完全反向

    # 斜向入射到yz平面
    ray_diagonal = Ray3D(PVector(0.0, 0.0, 0.0), PVector(1.0, 1.0, 0.0))
    plane_yz = Plane(PVector(1.0, 1.0, 0.0), PVector(1.0, -1.0, 0.0), PVector(1.0, 0.0, 1.0))  # x=1平面
    hit, point = PhysicalMeshes.intersect(ray_diagonal, plane_yz)
    @test hit
    @test point == PVector(1.0, 1.0, 0.0)
    ray_reflected_diag = reflect(ray_diagonal, plane_yz)
    @test ray_reflected_diag.n == PVector(-1.0, 1.0, 0.0)  # x方向反转
end

@testset "Ray2D - Boundary with units" begin
    # 带单位的平行测试
    ray_parallel_u = Ray2D(PVector2D(0.0u"m", 0.0u"m"), PVector2D(1.0u"m", 0.0u"m"))
    line_horizontal_u = Line2D(PVector2D(-1.0u"m", 1.0u"m"), PVector2D(1.0u"m", 1.0u"m"))
    hit, intersection = PhysicalMeshes.intersect(ray_parallel_u, line_horizontal_u)
    @test !hit
    @test isnothing(reflect(ray_parallel_u, line_horizontal_u))

    # 带单位的 miss 测试
    ray_right_u = Ray2D(PVector2D(0.0u"m", 0.0u"m"), PVector2D(1.0u"m", 0.0u"m"))
    line_left_u = Line2D(PVector2D(-2.0u"m", 1.0u"m"), PVector2D(-1.0u"m", 2.0u"m"))
    hit, intersection = PhysicalMeshes.intersect(ray_right_u, line_left_u)
    @test !hit
    @test isnothing(reflect(ray_right_u, line_left_u))

    # 带单位的45度反射测试
    ray_45_u = Ray2D(PVector2D(0.0u"m", 0.0u"m"), PVector2D(1.0u"m", 1.0u"m"))
    line_v_u = Line2D(PVector2D(2.0u"m", -1.0u"m"), PVector2D(2.0u"m", 3.0u"m"))
    hit, intersection = PhysicalMeshes.intersect(ray_45_u, line_v_u)
    @test hit
    @test intersection == PVector(2.0u"m", 2.0u"m")
    ray_reflected_u = reflect(ray_45_u, line_v_u)
    @test ray_reflected_u.theta ≈ 3π/4
end

@testset "Ray3D - Boundary with units" begin
    # 带单位的平行测试
    ray_parallel_u = Ray3D(PVector(0.0u"m", 0.0u"m", 0.0u"m"), PVector(1.0u"m", 0.0u"m", 0.0u"m"))
    plane_parallel_u = Plane(PVector(0.0u"m", 0.0u"m", 1.0u"m"), PVector(1.0u"m", 0.0u"m", 1.0u"m"), PVector(0.0u"m", 1.0u"m", 1.0u"m"))
    hit, point = PhysicalMeshes.intersect(ray_parallel_u, plane_parallel_u)
    @test !hit
    @test isnothing(reflect(ray_parallel_u, plane_parallel_u))

    # 带单位的 miss 测试
    ray_away_u = Ray3D(PVector(0.0u"m", 0.0u"m", 0.0u"m"), PVector(0.0u"m", 0.0u"m", -1.0u"m"))
    plane_above_u = Plane(PVector(1.0u"m", 1.0u"m", 1.0u"m"), PVector(1.0u"m", -1.0u"m", 1.0u"m"), PVector(-1.0u"m", 1.0u"m", 1.0u"m"))
    hit, point = PhysicalMeshes.intersect(ray_away_u, plane_above_u)
    @test !hit
    @test isnothing(reflect(ray_away_u, plane_above_u))

    # 带单位的45度反射测试
    ray_45_u = Ray3D(PVector(0.0u"m", 0.0u"m", 0.0u"m"), PVector(0.0u"m", 1.0u"m", 1.0u"m"))
    plane_xy_u = Plane(PVector(1.0u"m", 1.0u"m", 1.0u"m"), PVector(1.0u"m", -1.0u"m", 1.0u"m"), PVector(-1.0u"m", 1.0u"m", 1.0u"m"))
    hit, point = PhysicalMeshes.intersect(ray_45_u, plane_xy_u)
    @test hit
    @test point == PVector(0.0u"m", 1.0u"m", 1.0u"m")
    ray_reflected_u = reflect(ray_45_u, plane_xy_u)
    @test ray_reflected_u.n == PVector(0.0u"m", 1.0u"m", -1.0u"m")  # z方向反转
end

@testset "Normal vector validation" begin
    # 2D 垂直线法向量
    line_v = Line2D(PVector2D(2.0, 1.0), PVector2D(2.0, 3.0))  # 垂直线 x=2
    N_v = normal(line_v)
    @test N_v == PVector(-1.0, 0.0)  # 指向左侧

    # 2D 水平线法向量
    line_h = Line2D(PVector2D(1.0, 2.0), PVector2D(3.0, 2.0))  # 水平线 y=2
    N_h = normal(line_h)
    @test N_h == PVector(0.0, 1.0)  # 指向左侧法线方向

    # 2D 斜线法向量
    line_diag = Line2D(PVector2D(0.0, 0.0), PVector2D(1.0, 1.0))  # 对角线
    N_diag = normal(line_diag)
    @test N_diag ≈ PVector(-1.0, 1.0) / sqrt(2.0)  # 垂直于对角线 (单位化)

    # 3D 水平平面法向量
    plane_h = Plane(PVector(1.0, 1.0, 1.0), PVector(1.0, 1.0, -1.0), PVector(-1.0, 1.0, 1.0))  # y=1平面
    N_plane_h = normal(plane_h)
    @test N_plane_h == PVector(0.0, 1.0, 0.0)  # y方向

    # 3D 垂直平面法向量
    plane_v = Plane(PVector(1.0, 0.0, 0.0), PVector(1.0, 1.0, 0.0), PVector(1.0, 0.0, 1.0))  # x=1平面
    N_plane_v = normal(plane_v)
    @test N_plane_v == PVector(1.0, 0.0, 0.0)  # x方向
end
