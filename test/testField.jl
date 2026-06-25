using Test
using PhysicalParticles
using PhysicalMeshes

@testset "ArrayScalarField" begin
    # 构造函数
    field = ArrayScalarField(Float64, (10, 10, 10))
    @test field isa ArrayScalarField{Float64}
    @test size(field) == (10, 10, 10)
    @test length(field) == 1000
    @test all(field.data .== 0.0)

    # 从数组构造
    data = rand(5, 5, 5)
    field_from_array = ArrayScalarField(data)
    @test field_from_array.data == data

    # getindex
    field2 = ArrayScalarField(Float64, (5, 5, 5))
    field2.data[1, 1, 1] = 3.14
    @test field2[1, 1, 1] == 3.14
    @test field2[CartesianIndex(1, 1, 1)] == 3.14

    # setindex!
    field2[2, 2, 2] = 2.71
    @test field2.data[2, 2, 2] == 2.71

    # 边界访问
    field2[end, end, end] = 1.41
    @test field2[5, 5, 5] == 1.41

    # 显示
    @test occursin("ArrayScalarField", repr(field2))
end

@testset "ArrayVectorField" begin
    # 构造函数
    field = ArrayVectorField(Float64, (10, 10, 10), 3)
    @test field isa ArrayVectorField{Float64}
    @test size(field) == (10, 10, 10)
    @test length(field) == 1000
    @test field.dim == 3
    @test size(field.data) == (10, 10, 10, 3)

    # 2D 向量场
    field_2d = ArrayVectorField(Float64, (10, 10), 2)
    @test field_2d.dim == 2
    @test size(field_2d.data) == (10, 10, 2)

    # getindex - 返回向量
    field2 = ArrayVectorField(Float64, (5, 5, 5), 3)
    field2.data[1, 1, 1, :] = [1.0, 2.0, 3.0]
    vec = field2[1, 1, 1]
    @test vec == [1.0, 2.0, 3.0]

    # setindex! - 设置向量
    field2[2, 2, 2] = [4.0, 5.0, 6.0]
    @test field2.data[2, 2, 2, :] == [4.0, 5.0, 6.0]

    # 边界访问
    field2[end, end, end] = [7.0, 8.0, 9.0]
    @test field2[5, 5, 5] == [7.0, 8.0, 9.0]

    # 显示
    @test occursin("ArrayVectorField", repr(field2))
    @test occursin("dim 3", repr(field2))
end

@testset "ArrayTensorField" begin
    # 构造函数
    field = ArrayTensorField(Float64, (10, 10, 10), 3)
    @test field isa ArrayTensorField{Float64}
    @test size(field) == (10, 10, 10)
    @test length(field) == 1000
    @test field.dim == 3
    @test size(field.data) == (10, 10, 10, 3, 3)

    # 2D 张量场
    field_2d = ArrayTensorField(Float64, (10, 10), 2)
    @test field_2d.dim == 2
    @test size(field_2d.data) == (10, 10, 2, 2)

    # getindex - 返回矩阵
    field2 = ArrayTensorField(Float64, (5, 5, 5), 3)
    tensor = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    field2.data[1, 1, 1, :, :] = tensor
    mat = field2[1, 1, 1]
    @test mat == tensor

    # setindex! - 设置矩阵
    tensor2 = [9.0 8.0 7.0; 6.0 5.0 4.0; 3.0 2.0 1.0]
    field2[2, 2, 2] = tensor2
    @test field2.data[2, 2, 2, :, :] == tensor2

    # 边界访问
    tensor3 = ones(3, 3)
    field2[end, end, end] = tensor3
    @test field2[5, 5, 5] == tensor3

    # 显示
    @test occursin("ArrayTensorField", repr(field2))
    @test occursin("3x3", repr(field2))
end

@testset "Field types integration" begin
    # 创建一个网格并操作各种场
    mesh = MeshCartesianStatic(; Nx=10, Ny=10, Nz=10)

    # 标量场 - 密度
    rho = mesh.rho
    rho.data .= rand(size(rho.data)...)
    @test sum(rho.data) > 0

    # 向量场 - 速度
    vel = ArrayVectorField(Float64, size(rho), 3)
    for i in 1:size(vel, 1), j in 1:size(vel, 2), k in 1:size(vel, 3)
        vel[i, j, k] = [rand(), rand(), rand()]
    end
    @test size(vel.data)[end] == 3

    # 张量场 - 应力
    stress = ArrayTensorField(Float64, size(rho), 3)
    for i in 1:size(stress, 1), j in 1:size(stress, 2), k in 1:size(stress, 3)
        stress[i, j, k] = rand(3, 3)
    end
    @test size(stress.data)[end-1:end] == (3, 3)
end

@testset "Field boundary and error handling" begin
    # 标量场越界访问（应该抛出错误）
    scalar = ArrayScalarField(Float64, (5, 5, 5))
    @test_throws BoundsError scalar[6, 1, 1]
    @test_throws BoundsError scalar[0, 1, 1]

    # 向量场维度不匹配
    vector = ArrayVectorField(Float64, (5, 5, 5), 3)
    @test_throws DimensionMismatch vector[1, 1, 1] = [1.0, 2.0]  # 维度不匹配

    # 张量场维度不匹配
    tensor = ArrayTensorField(Float64, (5, 5, 5), 3)
    @test_throws DimensionMismatch tensor[1, 1, 1] = rand(2, 2)  # 维度不匹配
end

# =============================================================================
# Regression tests for lastindex / firstindex forwarding
#
# These guard against the bug surfaced by PhysicalFFT.jl CI on 2026-06-25:
# `MethodError: no method matching lastindex(::ArrayScalarField{...})` at the
# `m.phi[2:end-1]` site. The fix lives in FieldTypes.jl and forwards the
# AbstractArray iteration protocol to the underlying `data` field. The tests
# below exercise BOTH `end`-based and `2:end-1`-based indexing on every
# concrete field type, in 1D / 2D / 3D, so a regression is caught immediately.
# =============================================================================

@testset "Field: lastindex / firstindex (regression)" begin
    @testset "ArrayScalarField: 1D lastindex / firstindex" begin
        f = ArrayScalarField(Float64, (7,))
        @test lastindex(f) == 7
        @test firstindex(f) == 1
        @test lastindex(f, 1) == 7
        @test firstindex(f, 1) == 1

        # `end` resolution
        f.data .= collect(1.0:7.0)
        @test f[end] == 7.0

        # `2:end-1` slice (the exact pattern from PhysicalFFT.jl:342)
        f[2:end-1] .= 0.0
        @test f.data[2:6] == zeros(5)
        # Endpoints are preserved
        @test f.data[1] == 1.0
        @test f.data[7] == 7.0
    end

    @testset "ArrayScalarField: 2D lastindex / firstindex" begin
        f = ArrayScalarField(Float64, (4, 5))
        @test lastindex(f) == 20
        @test lastindex(f, 1) == 4
        @test lastindex(f, 2) == 5
        @test firstindex(f, 1) == 1
        @test firstindex(f, 2) == 1

        # `end` resolution
        @test f[end] == f[4, 5] == f.data[end]

        # Multi-dim `2:end-1` (PhysicalFFT.jl:346 pattern)
        f[2:end-1, 2:end-1] .= 1.0
        @test all(f.data[2:3, 2:4] .== 1.0)
        @test f.data[1, 1] == 0.0
        @test f.data[4, 5] == 0.0
    end

    @testset "ArrayScalarField: 3D lastindex / firstindex" begin
        f = ArrayScalarField(Float64, (3, 4, 5))
        @test lastindex(f) == 60
        @test lastindex(f, 1) == 3
        @test lastindex(f, 2) == 4
        @test lastindex(f, 3) == 5

        # 3D `2:end-1` (PhysicalFFT.jl:350 pattern)
        f[2:end-1, 2:end-1, 2:end-1] .= 2.0
        @test all(f.data[2:2, 2:3, 2:4] .== 2.0)
    end

    @testset "ArrayVectorField: spatial-only lastindex / firstindex" begin
        # `f.data` has trailing dim = vector dim (3), but `f` exposes only spatial dims.
        f = ArrayVectorField(Float64, (4, 5), 3)
        @test size(f) == (4, 5)
        @test size(f.data) == (4, 5, 3)
        @test length(f) == 20
        # The spatial-element count, not the data length
        @test lastindex(f) == 20
        @test lastindex(f, 1) == 4
        @test lastindex(f, 2) == 5
        @test firstindex(f) == 1
        @test firstindex(f, 1) == 1
        @test firstindex(f, 2) == 1

        # `end` returns the last *spatial* element as a vector
        f.data[end, end, :] = [7.0, 8.0, 9.0]
        @test f[end, end] == [7.0, 8.0, 9.0]

        # Multi-dim `2:end-1`: write the same 3-vector to every interior cell.
        # We use an explicit loop because `f[range, range]` flattens to a
        # (Nx*Ny, dim) matrix — broadcasting `[1,2,3]` into that shape is
        # semantically ambiguous. The `f[i, j] = v` form exercises both the
        # getindex/setindex! round-trip and the `lastindex(f, d)` resolution.
        for i in 2:3, j in 2:4
            f[i, j] = [1.0, 2.0, 3.0]
        end
        @test all(f.data[2:3, 2:4, :] .== reshape([1.0, 2.0, 3.0], 1, 1, 3))
    end

    @testset "ArrayVectorField: 3D spatial lastindex" begin
        f = ArrayVectorField(Float64, (3, 4, 5), 3)
        @test lastindex(f) == 60            # 3*4*5 spatial elements
        @test lastindex(f.data) == 60*3    # 3 trailing vector components
        @test lastindex(f, 1) == 3
        @test lastindex(f, 2) == 4
        @test lastindex(f, 3) == 5
    end

    @testset "ArrayTensorField: spatial-only lastindex / firstindex" begin
        # `f.data` has two trailing dims (3, 3) for tensor components.
        f = ArrayTensorField(Float64, (4, 5), 3)
        @test size(f) == (4, 5)
        @test size(f.data) == (4, 5, 3, 3)
        @test length(f) == 20
        @test lastindex(f) == 20
        @test lastindex(f, 1) == 4
        @test lastindex(f, 2) == 5
        @test firstindex(f) == 1

        # `end` returns the last spatial tensor (3x3 matrix)
        T = reshape(collect(1.0:9.0), 3, 3)
        f.data[end, end, :, :] = T
        @test f[end, end] == T

        # Multi-dim `2:end-1`: zero every interior tensor via setindex!.
        # Same rationale as the VectorField test above — the slice flattens
        # to (Nx*Ny, dim, dim) which doesn't have a natural broadcast target.
        Z = zeros(3, 3)
        for i in 2:3, j in 2:4
            f[i, j] = Z
        end
        @test all(f.data[2:3, 2:4, :, :] .== 0.0)
    end

    @testset "ArrayTensorField: 3D spatial lastindex" begin
        f = ArrayTensorField(Float64, (3, 4, 5), 3)
        @test lastindex(f) == 60
        @test lastindex(f, 3) == 5
    end
end

@testset "Field: integration with AbstractMesh (end-to-end)" begin
    # The original PhysicalFFT.jl bug surfaced via MeshCartesianStatic.phi;
    # verify that the mesh field indexing also works after the fix.
    m = MeshCartesianStatic(; Nx=5, Ny=5, Nz=5)  # default VertexMode → fields have size (6,6,6)
    @test m.phi isa ArrayScalarField
    @test lastindex(m.phi) == 6*6*6
    @test lastindex(m.phi, 1) == 6

    # The same slice pattern that crashed PhysicalFFT.jl precompile
    m.phi[2:end-1, 2:end-1, 2:end-1] .= 0.0
    @test all(m.phi.data[2:5, 2:5, 2:5] .== 0.0)

    # Boundary endpoints untouched
    @test m.phi.data[1, 1, 1] == 0.0
    @test m.phi.data[6, 6, 6] == 0.0
end
