@testset "Line" begin
    @testset "Line2D" begin
        line = Line2D(PVector2D(), PVector2D(1.0, 0.0))
        @test len(line) == 1.0
        @test midpoint(line) == PVector2D(0.5, 0.0)

        m = Line2D(PVector2D(2.0, 2.0, u"m"), PVector2D(2.0, 2.0, u"m"))
        n = Line2D(PVector2D(1.0, 1.0, u"m"), PVector2D(1.0, 1.0, u"m"))

        @test m + n == Line(PVector2D(3.0, 3.0, u"m"), PVector2D(3.0, 3.0, u"m"))
        @test m - n == Line(PVector2D(1.0, 1.0, u"m"), PVector2D(1.0, 1.0, u"m"))
        @test n * 2 == m
        @test 2 * n == m

        @test n + PVector2D(1.0, 1.0, u"m") == m
        @test m - PVector2D(1.0, 1.0, u"m") == n
        @test m / 2 == n
    end

    @testset "Line" begin
        line = Line(PVector(u"m"), PVector(1.0, 0.0, 0.0, u"m"))
        @test len(line) == 1.0u"m"
        @test midpoint(line) == PVector(0.5, 0.0, 0.0, u"m")

        m = Line(PVector(2.0, 2.0, 2.0), PVector(2.0, 2.0, 2.0))
        n = Line(PVector(1.0, 1.0, 1.0), PVector(1.0, 1.0, 1.0))

        @test m + n == Line(PVector(3.0, 3.0, 3.0), PVector(3.0, 3.0, 3.0))
        @test m - n == Line(PVector(1.0, 1.0, 1.0), PVector(1.0, 1.0, 1.0))
        @test n * 2 == m
        @test 2 * n == m

        @test n + PVector(1.0, 1.0, 1.0) == m
        @test m - PVector(1.0, 1.0, 1.0) == n
        @test m / 2 == n
    end

    # 边界测试
    @testset "Zero-length line" begin
        # Line2D 零长度测试
        line_zero_2d = Line2D(PVector2D(1.0, 1.0), PVector2D(1.0, 1.0))
        @test len(line_zero_2d) == 0.0
        @test midpoint(line_zero_2d) == PVector2D(1.0, 1.0)

        # Line 3D 零长度测试
        line_zero_3d = Line(PVector(1.0, 2.0, 3.0), PVector(1.0, 2.0, 3.0))
        @test len(line_zero_3d) == 0.0
        @test midpoint(line_zero_3d) == PVector(1.0, 2.0, 3.0)

        # 带单位的零长度线段
        line_zero_unit = Line(PVector(1.0, 1.0, 1.0, u"m"), PVector(1.0, 1.0, 1.0, u"m"))
        @test len(line_zero_unit) == 0.0u"m"
        @test midpoint(line_zero_unit) == PVector(1.0, 1.0, 1.0, u"m")
    end

    @testset "Negative coordinates" begin
        # 负坐标 Line2D
        line_neg_2d = Line2D(PVector2D(-5.0, -3.0), PVector2D(-1.0, -1.0))
        @test len(line_neg_2d) ≈ sqrt(16.0 + 4.0)
        @test midpoint(line_neg_2d) == PVector2D(-3.0, -2.0)

        # 负坐标 Line 3D
        line_neg_3d = Line(PVector(-1.0, -2.0, -3.0), PVector(-4.0, -5.0, -6.0))
        @test len(line_neg_3d) ≈ sqrt(9.0 + 9.0 + 9.0)
        @test midpoint(line_neg_3d) == PVector(-2.5, -3.5, -4.5)

        # 混合正负坐标
        line_mixed_2d = Line2D(PVector2D(-1.0, 2.0), PVector2D(3.0, -4.0))
        @test len(line_mixed_2d) ≈ sqrt(16.0 + 36.0)

        line_mixed_3d = Line(PVector(-1.0, 2.0, -3.0), PVector(4.0, -5.0, 6.0))
        @test len(line_mixed_3d) ≈ sqrt(25.0 + 49.0 + 81.0)
    end

    @testset "Extreme values" begin
        # 极大值
        line_large = Line(PVector(1e10, 1e10, 1e10), PVector(2e10, 2e10, 2e10))
        @test len(line_large) ≈ sqrt(3e20)

        line_large_2d = Line2D(PVector2D(1e10, 1e10), PVector2D(2e10, 2e10))
        @test len(line_large_2d) ≈ sqrt(2e20)

        # 极小值
        line_small = Line(PVector(1e-10, 1e-10, 1e-10), PVector(2e-10, 2e-10, 2e-10))
        @test len(line_small) ≈ sqrt(3e-20)

        line_small_2d = Line2D(PVector2D(1e-10, 1e-10), PVector2D(2e-10, 2e-10))
        @test len(line_small_2d) ≈ sqrt(2e-20)

        # 极大与极小混合（可能导致数值问题）
        line_mixed_extreme = Line(PVector(0.0, 0.0, 0.0), PVector(1e10, 1e-10, 0.0))
        @test len(line_mixed_extreme) ≈ sqrt(1e20 + 1e-20)
    end

    @testset "Arithmetic properties" begin
        line1 = Line2D(PVector2D(1.0, 2.0), PVector2D(3.0, 4.0))
        line2 = Line2D(PVector2D(5.0, 6.0), PVector2D(7.0, 8.0))

        # 加法交换律
        @test line1 + line2 == line2 + line1

        # 数乘分配律
        @test 2 * (line1 + line2) == 2 * line1 + 2 * line2

        # 与 PVector 运算
        vec = PVector2D(1.0, 1.0)
        @test (line1 + vec) - vec == line1

        # 3D 版本
        line3d_1 = Line(PVector(1.0, 2.0, 3.0), PVector(4.0, 5.0, 6.0))
        line3d_2 = Line(PVector(7.0, 8.0, 9.0), PVector(10.0, 11.0, 12.0))

        @test line3d_1 + line3d_2 == line3d_2 + line3d_1
        @test 3 * (line3d_1 + line3d_2) == 3 * line3d_1 + 3 * line3d_2

        vec3d = PVector(1.0, 1.0, 1.0)
        @test (line3d_1 + vec3d) - vec3d == line3d_1
    end

    @testset "Unit handling" begin
        # 不同单位组合
        line_m = Line(PVector(0.0, 0.0, 0.0, u"m"), PVector(1.0, 0.0, 0.0, u"m"))
        @test len(line_m) == 1.0u"m"

        line_km = Line(PVector(0.0, 0.0, 0.0, u"km"), PVector(1.0, 0.0, 0.0, u"km"))
        @test len(line_km) == 1.0u"km"

        # 2D 带单位
        line_2d_m = Line2D(PVector2D(0.0, 0.0, u"m"), PVector2D(3.0, 4.0, u"m"))
        @test len(line_2d_m) == 5.0u"m"

        # 标量乘法保持单位
        line_scaled = line_m * 2
        @test len(line_scaled) == 2.0u"m"

        # 除法保持单位
        line_divided = line_km / 2
        @test len(line_divided) == 0.5u"km"
    end

    @testset "Normal vector" begin
        # 2D 法向量测试
        line_horizontal = Line2D(PVector2D(0.0, 0.0), PVector2D(1.0, 0.0))
        n = normal(line_horizontal)
        @test n ≈ PVector2D(0.0, 1.0)  # 水平线的法向量应该是垂直向上

        line_vertical = Line2D(PVector2D(0.0, 0.0), PVector2D(0.0, 1.0))
        n = normal(line_vertical)
        @test n ≈ PVector2D(-1.0, 0.0)  # 垂直线的法向量应该是水平向左

        line_diagonal = Line2D(PVector2D(0.0, 0.0), PVector2D(1.0, 1.0))
        n = normal(line_diagonal)
        # 对角线的法向量应该是 (-1/√2, 1/√2) 或归一化后的值
        @test abs(n.x) ≈ abs(n.y)
        @test n.x ≈ -n.y

        # 负坐标法向量
        line_neg = Line2D(PVector2D(-1.0, -1.0), PVector2D(-2.0, -2.0))
        n = normal(line_neg)
        @test abs(n.x) ≈ abs(n.y)
    end
end
