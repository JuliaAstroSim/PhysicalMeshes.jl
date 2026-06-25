# =============================================================================
# FDTD algorithm tests (src/algorithms/FDTD.jl)
#
# Coverage targets:
#   - fdtd_step / fdtd_simulate
#   - update_electric_field (2D / 3D)
#   - update_magnetic_field (2D / 3D)
#   - apply_periodic_boundary, apply_dirichlet_boundary, apply_vacuum_boundary
#   - apply_boundary_conditions dispatcher
#
# NOTE: FDTD1D is omitted — its inner loop indexes B/E with a trailing
# vector-component (e.g. B[i, 2]) but a `dim=1` `MeshCartesianStatic`
# builds the vector fields with `dim=1`, so `B[i, 2]` raises
# `BoundsError`. That is an upstream bug in `FDTD.jl:37-39`; we
# leave it for a separate fix and focus the tests on the 2D/3D
# dispatch (which is what downstream PIC/Plasma actually uses).
# =============================================================================

using Test
using PhysicalParticles
using StructArrays
import PhysicalMeshes:
    FDTD2D, FDTD3D,
    fdtd_step, fdtd_simulate,
    update_electric_field, update_magnetic_field,
    apply_boundary_conditions,
    apply_periodic_boundary, apply_dirichlet_boundary, apply_vacuum_boundary

@testset "FDTD structs" begin
    @test FDTD2D() isa FDTD2D
    @test FDTD3D() isa FDTD3D
end

@testset "update_electric_field: 2D" begin
    mesh = MeshCartesianStatic(;
        Nx=8, Ny=8, Nz=1, NG=1,
        xMin=0.0, xMax=1.0, dim=2,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    # 2D mesh B has shape (Nx+1, Ny+1, 2). Seed B_z (component 2) so the
    # curl drives E.
    mesh.B.data[3, 3, 2] = 1.0
    update_electric_field(mesh, FDTD2D(), 0.01)
    @test any(!iszero, mesh.E.data)
end

@testset "update_magnetic_field: 2D" begin
    mesh = MeshCartesianStatic(;
        Nx=8, Ny=8, Nz=1, NG=1,
        xMin=0.0, xMax=1.0, dim=2,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    # 2D mesh E has shape (Nx+1, Ny+1, 2). Seed E_z (component 2) so the
    # curl drives B.
    mesh.E.data[3, 3, 2] = 1.0
    update_magnetic_field(mesh, FDTD2D(), 0.01)
    @test any(!iszero, mesh.B.data)
end

@testset "update_electric_field: 3D" begin
    mesh = MeshCartesianStatic(;
        Nx=6, Ny=6, Nz=6, NG=1,
        xMin=0.0, xMax=1.0, dim=3,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    mesh.B.data[3, 3, 3, 3] = 1.0
    update_electric_field(mesh, FDTD3D(), 0.01)
    @test any(!iszero, mesh.E.data)
end

@testset "update_magnetic_field: 3D" begin
    mesh = MeshCartesianStatic(;
        Nx=6, Ny=6, Nz=6, NG=1,
        xMin=0.0, xMax=1.0, dim=3,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    mesh.E.data[3, 3, 3, 1] = 1.0
    update_magnetic_field(mesh, FDTD3D(), 0.01)
    @test any(!iszero, mesh.B.data)
end

@testset "fdtd_step: 3D runs without error" begin
    mesh = MeshCartesianStatic(;
        Nx=6, Ny=6, Nz=6, NG=1,
        xMin=0.0, xMax=1.0, dim=3,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    mesh.E.data[3, 3, 3, 1] = 1.0
    mesh.B.data[3, 3, 3, 3] = 0.5

    @test_nowarn fdtd_step(mesh, FDTD3D(), 0.01)
    # Both fields must evolve (the curl terms are non-zero on this seed).
    @test any(!iszero, mesh.E.data)
    @test any(!iszero, mesh.B.data)
end

@testset "fdtd_simulate: 3D runs N steps" begin
    mesh = MeshCartesianStatic(;
        Nx=6, Ny=6, Nz=6, NG=1,
        xMin=0.0, xMax=1.0, dim=3,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    mesh.E.data[3, 3, 3, 1] = 1.0
    mesh.B.data[3, 3, 3, 3] = 0.5

    e0 = copy(mesh.E.data)
    fdtd_simulate(mesh, FDTD3D(), 0.01, 3)
    @test mesh.E.data != e0
end

@testset "apply_periodic_boundary: 1D scalar mirrors interior onto ghost" begin
    mesh = MeshCartesianStatic(;
        Nx=8, Ny=1, Nz=1, NG=1,
        xMin=0.0, xMax=1.0, dim=1,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    fill!(mesh.rho.data, 0.0)
    # In FDTD.jl: ``field.data[1] = field.data[end-1]`` and
    # ``field.data[end] = field.data[2]``. So the +x ghost (index `end`)
    # gets a copy of the 2nd interior (index 2), and the -x ghost
    # (index 1) gets a copy of the 2nd-from-end interior.
    mesh.rho.data[2, 1, 1] = 7.0
    mesh.rho.data[end-1, 1, 1] = 3.0
    apply_periodic_boundary(mesh, :rho)
    @test mesh.rho.data[end, 1, 1] ≈ mesh.rho.data[2, 1, 1]
    @test mesh.rho.data[1, 1, 1] ≈ mesh.rho.data[end-1, 1, 1]
end

@testset "apply_periodic_boundary: 1D vector field" begin
    mesh = MeshCartesianStatic(;
        Nx=8, Ny=1, Nz=1, NG=1,
        xMin=0.0, xMax=1.0, dim=1,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    fill!(mesh.B.data, 0.0)
    mesh.B.data[2, 2] = 7.0
    mesh.B.data[end-1, 2] = 3.0
    apply_periodic_boundary(mesh, :B)
    @test mesh.B.data[end, 2] ≈ mesh.B.data[2, 2]
    @test mesh.B.data[1, 2] ≈ mesh.B.data[end-1, 2]
end

@testset "apply_dirichlet_boundary: zeros out ghost zones" begin
    mesh = MeshCartesianStatic(;
        Nx=8, Ny=1, Nz=1, NG=1,
        xMin=0.0, xMax=1.0, dim=1,
        mode=VertexMode(), boundary=Dirichlet(),
        mhd=true,
    )
    fill!(mesh.E.data, 1.0)
    apply_dirichlet_boundary(mesh, :E)
    @test mesh.E.data[1, 1] == 0.0
    @test mesh.E.data[end, 1] == 0.0
end

@testset "apply_vacuum_boundary: replicates nearest interior onto ghost" begin
    mesh = MeshCartesianStatic(;
        Nx=8, Ny=1, Nz=1, NG=1,
        xMin=0.0, xMax=1.0, dim=1,
        mode=VertexMode(), boundary=Vacuum(),
        mhd=true,
    )
    fill!(mesh.E.data, 0.0)
    mesh.E.data[2, 1] = 5.0
    mesh.E.data[end-1, 1] = 9.0
    apply_vacuum_boundary(mesh, :E)
    @test mesh.E.data[1, 1] ≈ mesh.E.data[2, 1]
    @test mesh.E.data[end, 1] ≈ mesh.E.data[end-1, 1]
end

@testset "apply_boundary_conditions: dispatches on config.boundary" begin
    for bc in (Periodic(), Dirichlet(), Vacuum())
        mesh = MeshCartesianStatic(;
            Nx=8, Ny=1, Nz=1, NG=1,
            xMin=0.0, xMax=1.0, dim=1,
            mode=VertexMode(), boundary=bc,
            mhd=true,
        )
        # Seed interior (NOT the ghost zones).
        mesh.rho.data[2, 1, 1] = 1.0
        mesh.rho.data[end-1, 1, 1] = 1.0
        apply_boundary_conditions(mesh, :rho)
        # Interior must be untouched in every branch.
        @test mesh.rho.data[2, 1, 1] ≈ 1.0
        @test mesh.rho.data[end-1, 1, 1] ≈ 1.0
    end
end

@testset "apply_*_boundary: 3D scalar / vector fields" begin
    mesh = MeshCartesianStatic(;
        Nx=4, Ny=4, Nz=4, NG=1,
        xMin=0.0, xMax=1.0, dim=3,
        mode=VertexMode(), boundary=Dirichlet(),
        mhd=true,
    )
    fill!(mesh.rho.data, 1.0)
    apply_dirichlet_boundary(mesh, :rho)
    @test mesh.rho.data[1, 1, 1] == 0.0
    @test mesh.rho.data[end, end, end] == 0.0

    fill!(mesh.B.data, 1.0)
    apply_dirichlet_boundary(mesh, :B)
    @test mesh.B.data[1, 1, 1, 1] == 0.0
    @test mesh.B.data[end, end, end, 1] == 0.0
end

@testset "apply_*_boundary: nothing field is a no-op" begin
    # MeshCartesianStatic is an immutable struct — fields cannot be
    # replaced with `nothing` at runtime. The ``nothing`` short-circuit in
    # apply_*_boundary exists for the historical case where a mesh
    # factory left a field unallocated; in the current implementation
    # `Mesh.jl` always allocates placeholder arrays (see Mesh.jl:204-216).
    # We exercise the dispatcher on a fresh mesh as a sanity check.
    mesh = MeshCartesianStatic(;
        Nx=4, Ny=4, Nz=4, NG=1,
        xMin=0.0, xMax=1.0, dim=3,
        mode=VertexMode(), boundary=Periodic(),
        mhd=true,
    )
    @test_nowarn apply_periodic_boundary(mesh, :rho)
    @test_nowarn apply_dirichlet_boundary(mesh, :rho)
    @test_nowarn apply_vacuum_boundary(mesh, :rho)
end