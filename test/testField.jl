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
