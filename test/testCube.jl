@testset "Cube" begin
    c = Cube(PVector2D(1.0, 1.0), PVector2D())
    @test PhysicalMeshes.area(c) == 1.0

    c = Cube(PVector(1.0, 1.0, 1.0, u"m"), PVector(u"m"))
    @test PhysicalMeshes.volume(c) == 1.0u"m^3"

    @test interior(c, PVector(0.5, 0.5, 0.5, u"m"))
    @test exterior(c, PVector(1.5, 1.5, 1.5, u"m"))
    @test !interior(c, PVector(u"m"))
    @test !exterior(c, PVector(u"m"))

    # ========== 1. 零体积/面积测试 ==========
    @testset "Zero volume/area" begin
        # 2D 零面积 Cube
        c_zero_2d = Cube2D(PVector2D(1.0, 1.0), PVector2D(1.0, 1.0))
        @test PhysicalMeshes.area(c_zero_2d) == 0.0

        # 3D 零体积 Cube
        c_zero_3d = Cube(PVector(1.0, 1.0, 1.0), PVector(1.0, 1.0, 1.0))
        @test PhysicalMeshes.volume(c_zero_3d) == 0.0

        # 部分维度为零的 3D Cube
        c_zero_xy = Cube(PVector(1.0, 1.0, 2.0), PVector(1.0, 1.0, 1.0))
        @test PhysicalMeshes.volume(c_zero_xy) == 0.0

        c_zero_xz = Cube(PVector(1.0, 2.0, 1.0), PVector(1.0, 1.0, 1.0))
        @test PhysicalMeshes.volume(c_zero_xz) == 0.0

        c_zero_yz = Cube(PVector(2.0, 1.0, 1.0), PVector(1.0, 1.0, 1.0))
        @test PhysicalMeshes.volume(c_zero_yz) == 0.0
    end

    # ========== 2. 负坐标测试 ==========
    @testset "Negative coordinates" begin
        # 2D 负坐标
        c_neg_2d = Cube2D(PVector2D(-1.0, -1.0), PVector2D(-2.0, -2.0))
        @test PhysicalMeshes.area(c_neg_2d) == 1.0
        @test interior(c_neg_2d, PVector2D(-1.5, -1.5))
        @test !interior(c_neg_2d, PVector2D(0.0, 0.0))
        @test exterior(c_neg_2d, PVector2D(0.0, 0.0))
        @test !exterior(c_neg_2d, PVector2D(-1.5, -1.5))

        # 3D 负坐标
        c_neg_3d = Cube(PVector(-1.0, -1.0, -1.0), PVector(-2.0, -2.0, -2.0))
        @test PhysicalMeshes.volume(c_neg_3d) == 1.0
        @test interior(c_neg_3d, PVector(-1.5, -1.5, -1.5))
        @test !interior(c_neg_3d, PVector(0.0, 0.0, 0.0))
        @test exterior(c_neg_3d, PVector(0.0, 0.0, 0.0))

        # 跨越正负坐标的 Cube
        c_mixed = Cube(PVector(1.0, 1.0, 1.0), PVector(-1.0, -1.0, -1.0))
        @test PhysicalMeshes.volume(c_mixed) == 8.0
        @test interior(c_mixed, PVector(0.0, 0.0, 0.0))
        @test interior(c_mixed, PVector(-0.5, 0.5, -0.5))
    end

    # ========== 3. interior/exterior 边界测试 ==========
    @testset "Interior/Exterior boundary" begin
        c = Cube(PVector(2.0, 2.0, 2.0), PVector(0.0, 0.0, 0.0))

        # 内部点
        @test interior(c, PVector(1.0, 1.0, 1.0))
        @test interior(c, PVector(0.1, 0.1, 0.1))
        @test interior(c, PVector(1.9, 1.9, 1.9))
        @test interior(c, PVector(0.001, 1.0, 1.999))

        # 恰好在面上（边界情况）- interior 应该返回 false
        @test !interior(c, PVector(0.0, 1.0, 1.0))
        @test !interior(c, PVector(2.0, 1.0, 1.0))
        @test !interior(c, PVector(1.0, 0.0, 1.0))
        @test !interior(c, PVector(1.0, 2.0, 1.0))
        @test !interior(c, PVector(1.0, 1.0, 0.0))
        @test !interior(c, PVector(1.0, 1.0, 2.0))

        # 恰好在边上
        @test !interior(c, PVector(0.0, 0.0, 1.0))
        @test !interior(c, PVector(2.0, 2.0, 1.0))
        @test !interior(c, PVector(0.0, 2.0, 1.0))

        # 恰好在顶点
        @test !interior(c, PVector(0.0, 0.0, 0.0))
        @test !interior(c, PVector(2.0, 2.0, 2.0))
        @test !interior(c, PVector(0.0, 2.0, 2.0))
        @test !interior(c, PVector(2.0, 0.0, 0.0))

        # 外部点 - exterior 应该返回 true
        @test exterior(c, PVector(3.0, 3.0, 3.0))
        @test exterior(c, PVector(-1.0, 1.0, 1.0))
        @test exterior(c, PVector(1.0, -1.0, 1.0))
        @test exterior(c, PVector(1.0, 1.0, -1.0))
        @test exterior(c, PVector(2.1, 1.0, 1.0))

        # 内部点 - exterior 应该返回 false
        @test !exterior(c, PVector(1.0, 1.0, 1.0))
        @test !exterior(c, PVector(0.5, 0.5, 0.5))

        # 边界上 - exterior 应该返回 false（边界不算外部）
        @test !exterior(c, PVector(0.0, 1.0, 1.0))
        @test !exterior(c, PVector(2.0, 1.0, 1.0))
        @test !exterior(c, PVector(0.0, 0.0, 0.0))
    end

    # ========== 4. 不同尺寸测试 ==========
    @testset "Different sizes" begin
        # 单位 Cube
        c_unit = Cube(PVector(1.0, 1.0, 1.0), PVector(0.0, 0.0, 0.0))
        @test PhysicalMeshes.volume(c_unit) == 1.0

        # 非均匀尺寸
        c_rect = Cube(PVector(2.0, 3.0, 4.0), PVector(0.0, 0.0, 0.0))
        @test PhysicalMeshes.volume(c_rect) == 2.0 * 3.0 * 4.0
        @test PhysicalMeshes.volume(c_rect) == 24.0

        # 大 Cube
        c_large = Cube(PVector(100.0, 100.0, 100.0), PVector(0.0, 0.0, 0.0))
        @test PhysicalMeshes.volume(c_large) == 1e6

        # 小 Cube
        c_small = Cube(PVector(0.001, 0.001, 0.001), PVector(0.0, 0.0, 0.0))
        @test PhysicalMeshes.volume(c_small) == 1e-9
    end

    # ========== 5. 2D Cube (正方形) 测试 ==========
    @testset "2D Cube (Square)" begin
        # 基本 2D 测试
        c_2d = Cube2D(PVector2D(1.0, 1.0), PVector2D(0.0, 0.0))
        @test PhysicalMeshes.area(c_2d) == 1.0
        @test interior(c_2d, PVector2D(0.5, 0.5))
        @test !interior(c_2d, PVector2D(1.5, 0.5))
        @test exterior(c_2d, PVector2D(1.5, 0.5))

        # 2D 边界测试 - 边上
        @test !interior(c_2d, PVector2D(0.0, 0.5))
        @test !interior(c_2d, PVector2D(0.5, 0.0))
        @test !interior(c_2d, PVector2D(1.0, 0.5))
        @test !interior(c_2d, PVector2D(0.5, 1.0))

        # 2D 边界测试 - 顶点
        @test !interior(c_2d, PVector2D(0.0, 0.0))
        @test !interior(c_2d, PVector2D(1.0, 1.0))
        @test !interior(c_2d, PVector2D(0.0, 1.0))
        @test !interior(c_2d, PVector2D(1.0, 0.0))

        # 2D exterior 边界测试
        @test !exterior(c_2d, PVector2D(0.0, 0.5))
        @test !exterior(c_2d, PVector2D(0.0, 0.0))
        @test exterior(c_2d, PVector2D(-0.1, 0.5))
        @test exterior(c_2d, PVector2D(0.5, -0.1))

        # 2D 不同尺寸
        c_2d_rect = Cube2D(PVector2D(2.0, 3.0), PVector2D(0.0, 0.0))
        @test PhysicalMeshes.area(c_2d_rect) == 6.0

        c_2d_large = Cube2D(PVector2D(100.0, 100.0), PVector2D(0.0, 0.0))
        @test PhysicalMeshes.area(c_2d_large) == 10000.0
    end

    # ========== 6. 带单位测试增强 ==========
    @testset "Unitful support" begin
        # 带单位的 Cube
        c_unitful = Cube(PVector(2.0, 2.0, 2.0, u"m"), PVector(0.0, 0.0, 0.0, u"m"))
        @test PhysicalMeshes.volume(c_unitful) == 8.0u"m^3"
        @test interior(c_unitful, PVector(1.0, 1.0, 1.0, u"m"))
        @test !interior(c_unitful, PVector(2.5, 2.5, 2.5, u"m"))
        @test exterior(c_unitful, PVector(2.5, 2.5, 2.5, u"m"))

        # 2D 带单位
        c_2d_unitful = Cube2D(PVector2D(1.0, 1.0, u"m"), PVector2D(0.0, 0.0, u"m"))
        @test PhysicalMeshes.area(c_2d_unitful) == 1.0u"m^2"
        @test interior(c_2d_unitful, PVector2D(0.5, 0.5, u"m"))

        # 负坐标带单位
        c_neg_unitful = Cube(PVector(-1.0, -1.0, -1.0, u"m"), PVector(-2.0, -2.0, -2.0, u"m"))
        @test PhysicalMeshes.volume(c_neg_unitful) == 1.0u"m^3"
        @test interior(c_neg_unitful, PVector(-1.5, -1.5, -1.5, u"m"))
    end
end
