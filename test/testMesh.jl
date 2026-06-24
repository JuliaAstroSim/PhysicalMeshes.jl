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

    @testset "zoom function" begin
        # 1D zoom - 无插值
        a1 = [1, 2, 3]
        z1 = zoom(a1, 2)
        @test length(z1) == 6
        @test z1 == [1, 1, 2, 2, 3, 3]

        # 1D zoom - 有插值
        z1_interp = zoom(a1, 2, interp=LinearInterpolation)
        @test length(z1_interp) == 6
        @test z1_interp[1] ≈ 1.0
        @test z1_interp[end] ≈ 3.0

        # 2D zoom - 无插值
        a2 = [1 2; 3 4]
        z2 = zoom(a2, 2)
        @test size(z2) == (4, 4)
        @test z2[1, 1] == 1
        @test z2[3, 3] == 4

        # 2D zoom - 有插值
        z2_interp = zoom(a2, 2, interp=LinearInterpolation)
        @test size(z2_interp) == (4, 4)
        @test z2_interp[1, 1] ≈ 1.0
        @test z2_interp[end, end] ≈ 4.0

        # 3D zoom - 无插值
        a3 = ones(Int, 2, 2, 2)
        a3[1, 1, 1] = 1
        a3[2, 2, 2] = 8
        z3 = zoom(a3, 2)
        @test size(z3) == (4, 4, 4)
        @test z3[1, 1, 1] == 1
        @test z3[4, 4, 4] == 8

        # 3D zoom - 有插值
        z3_interp = zoom(a3, 2, interp=LinearInterpolation)
        @test size(z3_interp) == (4, 4, 4)
    end

    @testset "shrink function" begin
        # 1D shrink
        a1 = [1, 2, 3, 4, 5, 6]
        s1 = shrink(a1, 2)
        @test length(s1) == 3
        @test s1 == [2, 4, 6]

        # 1D shrink with shift
        s1_shift = shrink(a1, 2, 1)
        @test s1_shift == [1, 3, 5]

        # 2D shrink
        a2 = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]
        s2 = shrink(a2, 2)
        @test size(s2) == (2, 2)
        @test s2 == [6 8; 14 16]

        # 3D shrink
        a3 = reshape(collect(1:27), 3, 3, 3)
        s3 = shrink(a3, 3)
        @test size(s3) == (1, 1, 1)
        @test s3[1, 1, 1] == 27
    end

    @testset "mass_center function" begin
        # 创建测试网格
        mesh = MeshCartesianStatic(; Nx=10, Ny=10, Nz=10, xMin=-1.0, xMax=1.0)

        # 均匀密度 - 质心应该在中心
        mesh.rho.data .= 1.0
        center_uniform = mass_center(mesh)
        @test center_uniform.x ≈ 0.0 atol=0.1
        @test center_uniform.y ≈ 0.0 atol=0.1
        @test center_uniform.z ≈ 0.0 atol=0.1

        # 非均匀密度 - 质心偏向高密度区域
        mesh.rho.data .= 0.0
        # 在 x>0 区域设置高密度
        for i in 6:size(mesh.rho.data, 1)
            mesh.rho.data[i, :, :] .= 1.0
        end
        center_offset = mass_center(mesh)
        @test center_offset.x > 0.0

        # 零密度 - 应该返回几何中心
        mesh.rho.data .= 0.0
        center_zero = mass_center(mesh)
        # 验证返回的是几何中心（大致在0附近）
        @test abs(center_zero.x) < 2.0
        @test abs(center_zero.y) < 2.0
        @test abs(center_zero.z) < 2.0
    end

    @testset "total_mass function" begin
        # 创建测试网格
        mesh = MeshCartesianStatic(; Nx=10, Ny=10, Nz=10, xMin=0.0, xMax=1.0)
        cell_volume = prod(mesh.config.Δ)

        # 均匀密度 = 1
        mesh.rho.data .= 1.0
        total = total_mass(mesh)
        expected = sum(mesh.rho.data) * cell_volume
        @test total ≈ expected

        # 零密度
        mesh.rho.data .= 0.0
        @test total_mass(mesh) ≈ 0.0

        # 非均匀密度
        mesh.rho.data .= rand(size(mesh.rho.data)...)
        total_nonuniform = total_mass(mesh)
        expected_nonuniform = sum(mesh.rho.data) * cell_volume
        @test total_nonuniform ≈ expected_nonuniform
    end

    @testset "is_inbound function" begin
        config = MeshConfig(; Nx=10, Ny=10, Nz=10, xMin=-1.0, xMax=1.0)

        # 内部点
        @test is_inbound(PVector(0.0, 0.0, 0.0), config)
        @test is_inbound(PVector(0.5, 0.5, 0.5), config)
        @test is_inbound(PVector(-0.5, -0.5, -0.5), config)

        # 边界点（根据实现可能返回 true 或 false）
        @test is_inbound(PVector(1.0, 0.0, 0.0), config)
        @test is_inbound(PVector(-1.0, 0.0, 0.0), config)

        # 外部点
        @test !is_inbound(PVector(2.0, 0.0, 0.0), config)
        @test !is_inbound(PVector(0.0, 2.0, 0.0), config)
        @test !is_inbound(PVector(0.0, 0.0, 2.0), config)
    end

    @testset "outbound_list function" begin
        # 创建带粒子的网格 - 显式指定边界,避免 enlarge 将外侧粒子纳入
        pos = [
            PVector(0.0, 0.0, 0.0),   # 内部
            PVector(0.5, 0.5, 0.5),   # 内部
            PVector(2.0, 0.0, 0.0),   # 外部（x 超出）
            PVector(0.0, 2.0, 0.0),   # 外部（y 超出）
            PVector(0.0, 0.0, 2.0),   # 外部（z 超出）
        ]
        data = StructArray([Ball() for _ in 1:5])
        assign_particles(data, :Pos, pos)
        assign_particles(data, :Mass, 1.0)

        mesh = MeshCartesianStatic(data; Nx=10, Ny=10, Nz=10,
            xMin=-1.0, xMax=1.0, yMin=-1.0, yMax=1.0, zMin=-1.0, zMax=1.0,
            NG=0, enlarge=1.0, cube=false)

        # 获取越界列表
        outbound = outbound_list(mesh)
        @test 3 in outbound  # 第3个粒子越界
        @test 4 in outbound  # 第4个粒子越界
        @test 5 in outbound  # 第5个粒子越界
        @test !(1 in outbound)  # 第1个粒子在内部
        @test !(2 in outbound)  # 第2个粒子在内部
    end

    @testset "Mesh utility functions integration" begin
        # 创建一个实际场景：粒子分布到网格，然后计算质心和总质量
        n_particles = 100
        pos = [PVector(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1) for _ in 1:n_particles]
        data = StructArray([Ball() for _ in 1:n_particles])
        assign_particles(data, :Pos, pos)
        assign_particles(data, :Mass, 1.0)

        mesh = MeshCartesianStatic(data; Nx=20, Ny=20, Nz=20)

        # 验证质量守恒
        total_particle_mass = sum(data.Mass)
        total_mesh_mass = total_mass(mesh)
        @test total_mesh_mass ≈ total_particle_mass rtol=0.1

        # 使用 zoom 放大密度场
        rho_zoomed = zoom(mesh.rho.data, 2)
        @test size(rho_zoomed) == size(mesh.rho.data) .* 2

        # 使用 shrink 缩小密度场
        rho_shrunk = shrink(mesh.rho.data, 2)
        @test size(rho_shrunk) == size(mesh.rho.data) .÷ 2
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

    @testset "MeshCartesianStatic Modes" begin
        # CellMode - 属性位于单元中心
        config_cell = MeshConfig(; mode=CellMode(), Nx=10, Ny=10, Nz=10)
        mesh_cell = MeshCartesianStatic(config_cell)
        @test mesh_cell.config.mode isa CellMode
        # CellMode 下，数据点数量等于单元数量
        @test length(mesh_cell.rho.data) == 10 * 10 * 10

        # VertexMode - 属性位于网格点
        config_vertex = MeshConfig(; mode=VertexMode(), Nx=10, Ny=10, Nz=10)
        mesh_vertex = MeshCartesianStatic(config_vertex)
        @test mesh_vertex.config.mode isa VertexMode
        # VertexMode 下，数据点数量等于 (Nx+1) * (Ny+1) * (Nz+1)
        @test length(mesh_vertex.rho.data) == 11 * 11 * 11
    end

    @testset "MeshCartesianStatic Assignments" begin
        # 创建测试粒子
        pos = [PVector(0.0, 0.0, 0.0), PVector(0.5, 0.5, 0.5)]
        data = StructArray([Ball() for _ in 1:2])
        assign_particles(data, :Pos, pos)
        assign_particles(data, :Mass, 1.0)

        # NGP - Nearest Grid Point
        mesh_ngp = MeshCartesianStatic(data; assignment=NGP(), Nx=5, Ny=5, Nz=5)
        @test mesh_ngp.config.assignment isa NGP

        # CIC - Cloud In Cell
        mesh_cic = MeshCartesianStatic(data; assignment=CIC(), Nx=5, Ny=5, Nz=5)
        @test mesh_cic.config.assignment isa CIC

        # TSC - Triangular Shaped Cloud
        mesh_tsc = MeshCartesianStatic(data; assignment=TSC(), Nx=5, Ny=5, Nz=5)
        @test mesh_tsc.config.assignment isa TSC
    end

    @testset "MeshCartesianStatic Boundaries" begin
        # Periodic
        mesh_periodic = MeshCartesianStatic(; boundary=Periodic(), Nx=10, Ny=10, Nz=10)
        @test mesh_periodic.config.boundary isa Periodic

        # Vacuum
        mesh_vacuum = MeshCartesianStatic(; boundary=Vacuum(), Nx=10, Ny=10, Nz=10)
        @test mesh_vacuum.config.boundary isa Vacuum

        # Dirichlet
        mesh_dirichlet = MeshCartesianStatic(; boundary=Dirichlet(), Nx=10, Ny=10, Nz=10)
        @test mesh_dirichlet.config.boundary isa Dirichlet
    end

    @testset "MeshCartesianStatic Mass Conservation" begin
        # 创建均匀分布的粒子
        n_particles = 100
        pos = [PVector(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1) for _ in 1:n_particles]
        data = StructArray([Ball() for _ in 1:n_particles])
        assign_particles(data, :Pos, pos)
        assign_particles(data, :Mass, 1.0)

        total_particle_mass = sum(data.Mass)

        # 使用 CIC assignment 方法验证质量守恒
        mesh = MeshCartesianStatic(data; assignment=CIC(), Nx=20, Ny=20, Nz=20)
        total_mesh_mass = sum(mesh.rho.data) * prod(mesh.config.Δ)

        # 质量应该守恒（允许小的数值误差）
        @test total_mesh_mass ≈ total_particle_mass rtol=0.1
    end

    @testset "MeshCartesianStatic Constructors" begin
        # 从无构造
        mesh1 = MeshCartesianStatic()
        @test mesh1 isa MeshCartesianStatic

        # 从 MeshConfig 构造
        config = MeshConfig(; Nx=10, Ny=10, Nz=10)
        mesh2 = MeshCartesianStatic(config, nothing)
        @test mesh2 isa MeshCartesianStatic
        @test mesh2.config.N == config.N

        # 从粒子构造
        pos = [PVector(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1) for _ in 1:100]
        data = StructArray([Ball() for _ in 1:100])
        assign_particles(data, :Pos, pos)
        mesh3 = MeshCartesianStatic(data)
        @test mesh3 isa MeshCartesianStatic
        @test !isnothing(mesh3.data)
    end

    @testset "MeshCartesianStatic 2D" begin
        # 2D 网格
        mesh_2d = MeshCartesianStatic(; dim=2, Nx=10, Ny=10)
        @test mesh_2d.config.dim == 2
        @test length(mesh_2d.config.N) == 2
        @test length(mesh_2d.config.Δ) == 2
    end

    @testset "MeshCartesianStatic MHD Fields" begin
        mesh_mhd = MeshCartesianStatic(; mhd=true, Nx=10, Ny=10, Nz=10)

        # 验证 MHD 字段存在
        @test !isnothing(mesh_mhd.B)
        @test !isnothing(mesh_mhd.E)
        @test !isnothing(mesh_mhd.rho_e)
        @test !isnothing(mesh_mhd.j)

        # 验证字段尺寸与 pos 尺寸一致 (VertexMode => Len+1, CellMode => Len)
        expected_dims = mesh_mhd.config.mode isa VertexMode ?
            tuple((mesh_mhd.config.Len .+ 1)...) : tuple(mesh_mhd.config.Len...)
        @test size(mesh_mhd.B.data)[1:3] == expected_dims
        @test size(mesh_mhd.E.data)[1:3] == expected_dims
        @test size(mesh_mhd.rho_e.data) == expected_dims
        @test size(mesh_mhd.j.data)[1:3] == expected_dims
    end
end
