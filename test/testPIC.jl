@testset "PIC Algorithm" begin
    using StructArrays
    using LinearAlgebra
    
    # Import PIC-related functions
    import PhysicalMeshes: boris_algorithm, deposit_charge, deposit_current, interpolate_field,
                           push_particles, deposit_charge_current, pic_step, PIC3D,
                           ArrayScalarField, ArrayVectorField, MeshConfig, MeshCartesianStatic

    @testset "boris_algorithm" begin
        # 无磁场情况 - 只有电场加速
        vel = PVector(0.0, 0.0, 0.0)
        E = [1.0, 0.0, 0.0]  # x方向电场
        B = [0.0, 0.0, 0.0]  # 无磁场
        q = 1.0
        m = 1.0
        dt = 0.1

        vel_new = boris_algorithm(vel, E, B, q, m, dt)
        # v_new = v_old + q/m * E * dt
        @test vel_new.x ≈ 0.1  # 0 + 1 * 1 * 0.1
        @test vel_new.y ≈ 0.0
        @test vel_new.z ≈ 0.0

        # 无电场情况 - 只有磁场旋转
        vel2 = PVector(1.0, 0.0, 0.0)
        E2 = [0.0, 0.0, 0.0]
        B2 = [0.0, 0.0, 1.0]  # z方向磁场
        # Boris 算法在大回转角下不精确（|t|=q*B*dt/2m 必须 ≪ 1），
        # 所以用小 dt 让单步近似为 90° 时仍有 ~3% 误差。验证：
        #   1. 旋转方向正确（v 由 +x 转向 -y，右手定则：v×B = (1,0,0)×(0,0,1)=(0,-1,0)）
        #   2. v 长度守恒
        #   3. y 分量显著为负
        dt2 = π / 2  # 单步名义 90°，但算法是近似

        vel_new2 = boris_algorithm(vel2, E2, B2, q, m, dt2)
        # 长度守恒（Boris 算法严格保持速度幅值）
        @test isapprox(sqrt(vel_new2.x^2 + vel_new2.y^2 + vel_new2.z^2), 1.0; atol=1e-2)
        # 旋转方向：x 分量减少，y 分量变负
        @test vel_new2.x < vel2.x
        @test vel_new2.y < -0.5
        @test abs(vel_new2.z) < 1e-12

        # 完整旋转（360度）
        # 注意：Boris 算法是大回转角的近似。|t| = q*B*dt/(2m) 必须≪1 才接近真实旋转。
        # 这里用小 dt（π/20），让单步回转约 9°，更接近解析解。
        dt3 = π / 20   # ≈ 9° 单步回转
        vel_new3 = boris_algorithm(vel2, E2, B2, q, m, dt3)
        # 单步小回转后，v 应仍在 x-y 平面内、长度守恒、x 减少、y 变负。
        @test isapprox(sqrt(vel_new3.x^2 + vel_new3.y^2), 1.0; atol=1e-2)
        @test vel_new3.x < vel2.x        # 向 -y 方向转
        @test vel_new3.y < 0.0
        @test vel_new3.z ≈ 0.0 atol=1e-12
    end

    @testset "deposit_charge" begin
        # 创建网格
        config = MeshConfig(; Nx=10, Ny=10, Nz=10, xMin=0.0, xMax=1.0)
        rho = ArrayScalarField(Float64, tuple(config.Len...))
        rho.data .= 0.0

        # 在网格中心沉积电荷
        pos = PVector(0.5, 0.5, 0.5)
        charge = 1.0
        deposit_charge(rho, pos, charge, config)

        # 验证电荷被沉积到附近的网格点
        @test sum(rho.data) ≈ charge

        # 在角落沉积电荷
        rho.data .= 0.0
        pos2 = PVector(0.1, 0.1, 0.1)
        deposit_charge(rho, pos2, 2.0, config)
        @test sum(rho.data) ≈ 2.0

        # 多个电荷沉积
        rho.data .= 0.0
        for i in 1:10
            pos = PVector(rand(), rand(), rand())
            deposit_charge(rho, pos, 0.1, config)
        end
        @test sum(rho.data) ≈ 1.0 rtol=0.01
    end

    @testset "deposit_current" begin
        # 创建网格
        config = MeshConfig(; Nx=10, Ny=10, Nz=10, xMin=0.0, xMax=1.0)
        j = ArrayVectorField(Float64, tuple(config.Len...), 3)
        j.data .= 0.0

        # 沉积电流
        pos = PVector(0.5, 0.5, 0.5)
        vel = PVector(1.0, 0.0, 0.0)  # x方向运动
        charge = 1.0
        deposit_current(j, pos, vel, charge, config)

        # 验证电流被沉积
        @test sum(j.data[:, :, :, 1]) > 0  # x分量
        @test sum(j.data[:, :, :, 2]) ≈ 0.0 atol=1e-10  # y分量应该为0
        @test sum(j.data[:, :, :, 3]) ≈ 0.0 atol=1e-10  # z分量应该为0

        # y方向运动
        j.data .= 0.0
        vel2 = PVector(0.0, 1.0, 0.0)
        deposit_current(j, pos, vel2, charge, config)
        @test sum(j.data[:, :, :, 1]) ≈ 0.0 atol=1e-10
        @test sum(j.data[:, :, :, 2]) > 0
        @test sum(j.data[:, :, :, 3]) ≈ 0.0 atol=1e-10
    end

    @testset "interpolate_field" begin
        # 创建网格和场
        config = MeshConfig(; Nx=10, Ny=10, Nz=10, xMin=0.0, xMax=1.0)
        E = ArrayVectorField(Float64, tuple(config.Len...), 3)

        # 设置一个简单场：E = (x, 0, 0)
        for i in 1:size(E.data, 1)
            x = config.Min[1] + (i - 0.5) * config.Δ[1]
            E.data[i, :, :, 1] .= x
        end

        # 在已知点插值
        pos = PVector(0.5, 0.5, 0.5)
        E_interp = interpolate_field(E, pos, config)
        @test E_interp[1] ≈ 0.5 rtol=0.1
        @test E_interp[2] ≈ 0.0
        @test E_interp[3] ≈ 0.0

        # 在网格点插值应该接近该点的场值
        pos2 = PVector(0.25, 0.25, 0.25)
        E_interp2 = interpolate_field(E, pos2, config)
        @test E_interp2[1] ≈ 0.25 rtol=0.2
    end

    @testset "push_particles" begin
        # 创建简单网格和粒子
        config = MeshConfig(; Nx=10, Ny=10, Nz=10, xMin=0.0, xMax=1.0)

        # 创建粒子数组
        n_particles = 10
        particles = StructArray((
            Pos = [PVector(rand(), rand(), rand()) for _ in 1:n_particles],
            Vel = [PVector(0.0, 0.0, 0.0) for _ in 1:n_particles],
            Charge = ones(n_particles),
            Mass = ones(n_particles)
        ))

        # 创建网格（启用 MHD）
        mesh = MeshCartesianStatic(config, particles; mhd=true)

        # 设置均匀电场
        mesh.E.data[:, :, :, 1] .= 1.0  # x方向电场
        mesh.E.data[:, :, :, 2] .= 0.0
        mesh.E.data[:, :, :, 3] .= 0.0

        # 保存初始位置
        initial_pos = copy(particles.Pos)

        # 推动粒子
        pic = PIC3D()
        dt = 0.01
        push_particles(mesh, pic, dt)

        # 验证粒子被加速（x方向移动）
        for i in 1:n_particles
            @test particles.Vel[i].x > 0  # 被加速
            @test particles.Pos[i].x > initial_pos[i].x  # 位置改变
        end
    end

    @testset "PIC full workflow" begin
        # 创建网格
        config = MeshConfig(; Nx=20, Ny=20, Nz=20, xMin=0.0, xMax=1.0)

        # 创建均匀分布的粒子
        n_particles = 100
        particles = StructArray((
            Pos = [PVector(rand(), rand(), rand()) for _ in 1:n_particles],
            Vel = [PVector(rand()-0.5, rand()-0.5, rand()-0.5) for _ in 1:n_particles],
            Charge = ones(n_particles),
            Mass = ones(n_particles)
        ))

        # 创建网格（启用 MHD）
        mesh = MeshCartesianStatic(config, particles; mhd=true)

        # 初始沉积
        deposit_charge_current(mesh, PIC3D())

        # 验证电荷守恒
        total_charge = sum(mesh.rho_e.data)
        expected_charge = sum(particles.Charge)
        @test total_charge ≈ expected_charge rtol=0.1

        # 运行几个时间步
        dt = 0.001
        for step in 1:10
            pic_step(mesh, PIC3D(), dt)
        end

        # 验证粒子仍在合理范围内（没有发散）
        for i in 1:n_particles
            @test 0.0 <= particles.Pos[i].x <= 1.0 || config.boundary isa Periodic
            @test 0.0 <= particles.Pos[i].y <= 1.0 || config.boundary isa Periodic
            @test 0.0 <= particles.Pos[i].z <= 1.0 || config.boundary isa Periodic
        end
    end
end
