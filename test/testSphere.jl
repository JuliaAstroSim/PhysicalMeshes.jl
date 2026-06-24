@testset "Sphere" begin
    s = Sphere(PVector(), 1.0)
    @test interior(s, PVector(0.1, 0.1, 0.1))
    @test exterior(s, PVector(1.0, 1.0, 1.0))
    @test !interior(s, PVector(1.0, 0.0, 0.0))
    @test !exterior(s, PVector(1.0, 0.0, 0.0))

    # 1. 零半径测试
    @testset "Zero radius" begin
        s_zero = Sphere(PVector(0.0, 0.0, 0.0), 0.0)
        # 零半径球体，只有中心点可能在内部（距离 = 0 < 0 为 false）
        @test !interior(s_zero, PVector(0.0, 0.0, 0.0))
        # 中心点也不是外部（距离 = 0 > 0 为 false）
        @test !exterior(s_zero, PVector(0.0, 0.0, 0.0))
        # 其他点都是外部
        @test exterior(s_zero, PVector(0.001, 0.0, 0.0))
        @test exterior(s_zero, PVector(1.0, 0.0, 0.0))
    end

    # 2. interior/exterior 边界测试
    @testset "Interior/Exterior boundary" begin
        s = Sphere(PVector(0.0, 0.0, 0.0), 1.0)

        # 内部点
        @test interior(s, PVector(0.0, 0.0, 0.0))  # 中心
        @test interior(s, PVector(0.5, 0.0, 0.0))
        @test interior(s, PVector(0.0, 0.5, 0.0))
        @test interior(s, PVector(0.0, 0.0, 0.5))
        @test interior(s, PVector(0.5, 0.5, 0.5))
        @test interior(s, PVector(-0.5, 0.0, 0.0))
        @test interior(s, PVector(0.0, -0.5, 0.0))
        @test interior(s, PVector(0.0, 0.0, -0.5))

        # 恰好在球面上（边界情况）- 既不是 interior 也不是 exterior
        @test !interior(s, PVector(1.0, 0.0, 0.0))
        @test !exterior(s, PVector(1.0, 0.0, 0.0))
        @test !interior(s, PVector(0.0, 1.0, 0.0))
        @test !exterior(s, PVector(0.0, 1.0, 0.0))
        @test !interior(s, PVector(0.0, 0.0, 1.0))
        @test !exterior(s, PVector(0.0, 0.0, 1.0))
        @test !interior(s, PVector(-1.0, 0.0, 0.0))
        @test !exterior(s, PVector(-1.0, 0.0, 0.0))

        # 球面上的点（对角线方向）
        r = 1.0 / sqrt(3.0)
        @test !interior(s, PVector(r, r, r))
        @test !exterior(s, PVector(r, r, r))

        # 外部点
        @test exterior(s, PVector(2.0, 0.0, 0.0))
        @test exterior(s, PVector(0.0, 2.0, 0.0))
        @test exterior(s, PVector(0.0, 0.0, 2.0))
        @test exterior(s, PVector(1.0, 1.0, 1.0))
        @test exterior(s, PVector(-2.0, 0.0, 0.0))

        # 非外部点（内部或边界）
        @test !exterior(s, PVector(0.0, 0.0, 0.0))
        @test !exterior(s, PVector(0.5, 0.0, 0.0))
    end

    # 3. 不同半径测试
    @testset "Different radii" begin
        # 小半径
        s_small = Sphere(PVector(0.0, 0.0, 0.0), 0.001)
        @test interior(s_small, PVector(0.0, 0.0, 0.0))
        @test exterior(s_small, PVector(0.01, 0.0, 0.0))
        @test !interior(s_small, PVector(0.001, 0.0, 0.0))  # 边界
        @test !exterior(s_small, PVector(0.001, 0.0, 0.0))  # 边界

        # 大半径
        s_large = Sphere(PVector(0.0, 0.0, 0.0), 1e6)
        @test interior(s_large, PVector(1e5, 1e5, 1e5))
        @test exterior(s_large, PVector(2e6, 0.0, 0.0))
        @test interior(s_large, PVector(5e5, 5e5, 5e5))

        # 中等半径
        s_medium = Sphere(PVector(0.0, 0.0, 0.0), 10.0)
        @test interior(s_medium, PVector(5.0, 5.0, 5.0))
        @test !interior(s_medium, PVector(10.0, 0.0, 0.0))  # 边界
        @test exterior(s_medium, PVector(15.0, 0.0, 0.0))
    end

    # 4. 不同中心位置测试
    @testset "Different center positions" begin
        # 中心不在原点
        s_offset = Sphere(PVector(5.0, 5.0, 5.0), 1.0)
        @test interior(s_offset, PVector(5.0, 5.0, 5.0))  # 中心
        @test interior(s_offset, PVector(5.5, 5.0, 5.0))
        @test interior(s_offset, PVector(5.0, 5.5, 5.0))
        @test interior(s_offset, PVector(5.0, 5.0, 5.5))
        @test exterior(s_offset, PVector(0.0, 0.0, 0.0))
        @test exterior(s_offset, PVector(7.0, 5.0, 5.0))
        @test exterior(s_offset, PVector(5.0, 7.0, 5.0))

        # 负坐标中心
        s_neg = Sphere(PVector(-5.0, -5.0, -5.0), 2.0)
        @test interior(s_neg, PVector(-5.0, -5.0, -5.0))
        @test interior(s_neg, PVector(-4.0, -5.0, -5.0))
        @test interior(s_neg, PVector(-5.0, -4.0, -5.0))
        @test exterior(s_neg, PVector(0.0, 0.0, 0.0))
        @test exterior(s_neg, PVector(-8.0, -5.0, -5.0))

        # 混合坐标中心
        s_mixed = Sphere(PVector(-3.0, 4.0, -2.0), 3.0)
        @test interior(s_mixed, PVector(-3.0, 4.0, -2.0))
        @test interior(s_mixed, PVector(-2.0, 4.0, -2.0))
        @test exterior(s_mixed, PVector(0.0, 0.0, 0.0))
    end

    # 5. 带单位测试增强
    @testset "Unitful spheres" begin
        # 带单位的 Sphere
        s_unitful = Sphere(PVector(0.0, 0.0, 0.0, u"m"), 1.0u"m")
        @test interior(s_unitful, PVector(0.0, 0.0, 0.0, u"m"))
        @test interior(s_unitful, PVector(0.5, 0.0, 0.0, u"m"))
        @test exterior(s_unitful, PVector(2.0, 0.0, 0.0, u"m"))
        @test !interior(s_unitful, PVector(1.0, 0.0, 0.0, u"m"))  # 边界

        # 不同单位
        s_km = Sphere(PVector(0.0, 0.0, 0.0, u"km"), 1.0u"km")
        @test interior(s_km, PVector(0.5, 0.0, 0.0, u"km"))
        @test exterior(s_km, PVector(1.5, 0.0, 0.0, u"km"))

        # 天文单位
        s_au = Sphere(PVector(0.0, 0.0, 0.0, u"AU"), 1.0u"AU")
        @test interior(s_au, PVector(0.5, 0.0, 0.0, u"AU"))
        @test exterior(s_au, PVector(1.5, 0.0, 0.0, u"AU"))
    end

    # 6. 边界一致性测试
    @testset "Boundary consistency" begin
        s = Sphere(PVector(0.0, 0.0, 0.0), 1.0)

        # 验证 interior 和 exterior 互斥（除了边界）
        test_points = [
            (-0.9, 0.0, 0.0), (-0.5, 0.0, 0.0), (0.0, 0.0, 0.0),
            (0.5, 0.0, 0.0), (0.9, 0.0, 0.0),  # 内部
            (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),  # 边界
            (1.1, 0.0, 0.0), (2.0, 0.0, 0.0),  # 外部
            (0.5, 0.5, 0.0), (0.5, 0.5, 0.5),
            (-0.5, -0.5, 0.0), (-0.5, -0.5, -0.5),
        ]

        for (x, y, z) in test_points
            p = PVector(x, y, z)
            r = sqrt(x^2 + y^2 + z^2)
            if r < 1.0
                @test interior(s, p)
                @test !exterior(s, p)
            elseif r > 1.0
                @test !interior(s, p)
                @test exterior(s, p)
            else
                # 边界情况：既不是 interior 也不是 exterior
                @test !interior(s, p)
                @test !exterior(s, p)
            end
        end
    end
end
