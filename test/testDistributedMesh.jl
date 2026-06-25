# =============================================================================
# DistributedMesh tests (src/parallel/DistributedMesh.jl)
#
# Coverage targets:
#   - partition_neighbors (1D / 2D / 3D)
#   - calculate_local_config (1D / 2D / 3D)
#   - Base.getproperty forwarding to local mesh
#   - DistributedMesh(config, units) constructor
#   - DistributedMesh(config, local_mesh, rank, size, neighbors) explicit constructor
#   - Base.show
#   - AbstractMesh fallbacks (number_of_nodes / cells / fields)
# =============================================================================

using Test
using PhysicalParticles
using StaticArrays
import PhysicalMeshes:
    DistributedMesh,
    partition_neighbors, calculate_local_config,
    AbstractMesh, MeshConfig

@testset "partition_neighbors: 1D chain" begin
    # Linear layout: rank 0 only has neighbour 1; rank in middle has both
    # neighbours; last rank only has its left neighbour.
    @test partition_neighbors(0, 2, 1) == [1]
    @test partition_neighbors(1, 2, 1) == [0]
    @test partition_neighbors(1, 3, 1) == [0, 2]
    @test partition_neighbors(0, 4, 1) == [1]
    @test partition_neighbors(3, 4, 1) == [2]
    @test partition_neighbors(2, 5, 1) == [1, 3]
end

@testset "partition_neighbors: 2D grid (sqrt layout)" begin
    # 2D uses rows = ceil(sqrt(size)), cols = ceil(size / rows).
    # 4 ranks in 2x2 grid: every interior rank has 4 neighbours; corner ranks
    # have 2 neighbours.
    @test sort(partition_neighbors(0, 4, 2)) == [1, 2]
    @test sort(partition_neighbors(1, 4, 2)) == [0, 3]
    @test sort(partition_neighbors(2, 4, 2)) == [0, 3]
    @test sort(partition_neighbors(3, 4, 2)) == [1, 2]
end

@testset "partition_neighbors: 3D (uses 1D layout)" begin
    # 3D partition uses the same linear chain as 1D: neighbours are
    # rank-1 and rank+1 when they exist.
    @test partition_neighbors(0, 3, 3) == [1]
    @test partition_neighbors(1, 3, 3) == [0, 2]
    @test partition_neighbors(2, 3, 3) == [1]
end

@testset "partition_neighbors: every returned rank is in [0, size)" begin
    for dim in (1, 2, 3), size in 1:6
        for rank in 0:(size-1)
            for n in partition_neighbors(rank, size, dim)
                @test 0 <= n < size
            end
        end
    end
end

@testset "calculate_local_config: 1D divides Nx across ranks" begin
    base = MeshConfig(;
        Nx=10, Ny=1, Nz=1, dim=1,
        xMin=0.0, xMax=1.0,
        mode=VertexMode(), boundary=Periodic(),
    )

    # 1 rank → that rank owns the whole domain.
    cfg0 = calculate_local_config(base, 0, 1)
    @test cfg0.N == SVector(10)
    @test cfg0.Min[1] ≈ 0.0
    @test cfg0.Max[1] ≈ 1.0

    # 2 ranks → split Nx=10 into 5+5 (no remainder).
    cfg0_2 = calculate_local_config(base, 0, 2)
    cfg1_2 = calculate_local_config(base, 1, 2)
    @test cfg0_2.N[1] + cfg1_2.N[1] == 10
    @test cfg0_2.Max[1] ≈ cfg1_2.Min[1]   # adjacency

    # 3 ranks → 10/3 = 3, remainder = 1, so ranks 0,1,2 get 4,3,3 cells.
    cfgs = [calculate_local_config(base, r, 3) for r in 0:2]
    @test [cfgs[r+1].N[1] for r in 0:2] == [4, 3, 3]
    @test sum(c.N[1] for c in cfgs) == 10
end

@testset "calculate_local_config: 2D row/column layout" begin
    base = MeshConfig(;
        Nx=4, Ny=4, Nz=1, dim=2,
        xMin=0.0, xMax=1.0, yMin=0.0, yMax=1.0,
        mode=VertexMode(), boundary=Periodic(),
    )

    # 4 ranks → 2x2 grid → each rank owns 2x2 cells.
    for r in 0:3
        cfg = calculate_local_config(base, r, 4)
        @test cfg.N[1] == 2
        @test cfg.N[2] == 2
    end

    # Invariant: along each axis the local ranks tile the global extent.
    # We check that rank i+1 in a column picks up where rank i left off
    # (i.e., the local Min[1] of rank (row, col+1) == local Max[1] of
    # (row, col)).
    cfgs = [calculate_local_config(base, r, 4) for r in 0:3]
    # Row 0: ranks 0 and 1 → adjacency in x
    @test cfgs[1].Min[1] == 0.0
    @test cfgs[2].Min[1] ≈ cfgs[1].Max[1]
    # Column 0: ranks 0 and 2 → adjacency in y
    @test cfgs[3].Min[2] ≈ cfgs[1].Max[2]
end

@testset "calculate_local_config: 3D decomposes only along x" begin
    base = MeshConfig(;
        Nx=8, Ny=4, Nz=2, dim=3,
        xMin=0.0, xMax=1.0,
        mode=VertexMode(), boundary=Periodic(),
    )
    # 2 ranks → split Nx=8 into 4+4, Ny/Nz preserved.
    cfg0 = calculate_local_config(base, 0, 2)
    cfg1 = calculate_local_config(base, 1, 2)
    @test cfg0.N == SVector(4, 4, 2)
    @test cfg1.N == SVector(4, 4, 2)
    @test cfg0.Max[2] ≈ cfg1.Max[2]   # y is shared
    @test cfg0.Max[3] ≈ cfg1.Max[3]   # z is shared
end

@testset "calculate_local_config: preserves dim, mode, boundary" begin
    base = MeshConfig(;
        Nx=8, Ny=8, Nz=1, dim=2,
        mode=CellMode(), boundary=Dirichlet(),
        xMin=0.0, xMax=1.0, yMin=0.0, yMax=1.0,
    )
    cfg = calculate_local_config(base, 0, 4)
    @test cfg.dim == 2
    @test cfg.mode isa CellMode
    @test cfg.boundary isa Dirichlet
end

@testset "DistributedMesh: explicit constructor stores rank/size/neighbours" begin
    cfg = MeshConfig(;
        Nx=4, Ny=4, Nz=4, dim=3,
        xMin=0.0, xMax=1.0, boundary=Periodic(),
        mode=VertexMode(),
    )
    local_mesh = MeshCartesianStatic(cfg, nothing, nothing)
    dm = DistributedMesh(cfg, local_mesh, 0, 2, [1])

    @test dm isa DistributedMesh
    @test dm isa AbstractMesh
    @test dm.rank == 0
    @test dm.size == 2
    @test dm.neighbors == [1]
    @test dm.config === cfg
end

@testset "DistributedMesh: getproperty forwards to local_mesh" begin
    cfg = MeshConfig(;
        Nx=4, Ny=4, Nz=4, dim=3,
        xMin=0.0, xMax=1.0, boundary=Periodic(),
        mode=VertexMode(),
    )
    local_mesh = MeshCartesianStatic(cfg, nothing, nothing)
    local_mesh.rho.data .= 5.0
    dm = DistributedMesh(cfg, local_mesh, 0, 1, Int[])

    # Field access through DistributedMesh should land on the local mesh.
    @test dm.rho.data[1, 1, 1] ≈ 5.0
    @test dm.config === cfg
    @test dm.local_mesh === local_mesh
end

@testset "DistributedMesh: number_of_nodes/cells delegate to local_mesh" begin
    cfg = MeshConfig(;
        Nx=4, Ny=4, Nz=4, dim=3,
        xMin=0.0, xMax=1.0, boundary=Periodic(),
        mode=VertexMode(),
    )
    local_mesh = MeshCartesianStatic(cfg, nothing, nothing)
    dm = DistributedMesh(cfg, local_mesh, 0, 1, Int[])

    @test number_of_nodes(dm) == number_of_nodes(local_mesh)
    @test number_of_cells(dm) == number_of_cells(local_mesh)
    @test number_of_fields(dm) == number_of_fields(local_mesh)
end

@testset "DistributedMesh: Base.show does not error" begin
    cfg = MeshConfig(;
        Nx=4, Ny=4, Nz=4, dim=3,
        xMin=0.0, xMax=1.0, boundary=Periodic(),
        mode=VertexMode(),
    )
    local_mesh = MeshCartesianStatic(cfg, nothing, nothing)
    dm = DistributedMesh(cfg, local_mesh, 0, 2, [1])
    s = sprint(show, dm)
    @test occursin("DistributedMesh", s)
    @test occursin("rank", s)
end

@testset "DistributedMesh: convenience constructor (master process, no workers)" begin
    # The convenience constructor places the master process at rank 0
    # (`rank = myid() - 1`) and computes `size = nworkers() + 1`. In a
    # single-process CI run the master pid may or may not appear in
    # `workers()`, so we just check that the rank and size are sane and
    # the constructor returns a `DistributedMesh`.
    cfg = MeshConfig(;
        Nx=4, Ny=4, Nz=4, dim=3,
        xMin=0.0, xMax=1.0, boundary=Periodic(),
        mode=VertexMode(),
    )
    dm = DistributedMesh(cfg)
    @test dm isa DistributedMesh
    @test dm.rank >= 0
    @test dm.size >= 1
    # The local mesh must be a valid MeshCartesianStatic.
    @test dm.local_mesh isa MeshCartesianStatic
    # The master rank always owns the start of the domain.
    @test dm.local_mesh.config.Min[1] ≈ cfg.Min[1]
end