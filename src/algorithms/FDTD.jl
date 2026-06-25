abstract type AbstractFDTD end

struct FDTD1D <: AbstractFDTD end
struct FDTD2D <: AbstractFDTD end
struct FDTD3D <: AbstractFDTD end

function fdtd_step(mesh::AbstractMesh, fdtd::AbstractFDTD, dt::Float64)
    # Update electric field
    update_electric_field(mesh, fdtd, dt)

    # Synchronize electric field (for distributed meshes)
    if mesh isa DistributedMesh
        sync_field(mesh, :E)
    end

    # Update magnetic field
    update_magnetic_field(mesh, fdtd, dt)

    # Synchronize magnetic field (for distributed meshes)
    if mesh isa DistributedMesh
        sync_field(mesh, :B)
    end
end

function update_electric_field(mesh::AbstractMesh, fdtd::FDTD1D, dt::Float64)
    # 1D FDTD electric field update
    config = mesh.config
    E = mesh.E
    B = mesh.B
    j = mesh.j

    dx = config.Δ[1]
    cdt_dx = dt / dx

    # Update E_x only (1D)
    @inbounds for i in 2:size(E, 1)-1
        E[i, 1] += cdt_dx * (B[i+1, 2] - B[i, 2])  # dB_z/dx (z component is index 2 in 1D)
        E[i, 1] -= dt * j[i, 1]
    end

    apply_boundary_conditions(mesh, :E)
end

function update_electric_field(mesh::AbstractMesh, fdtd::FDTD2D, dt::Float64)
    # 2D FDTD electric field update
    config = mesh.config
    E = mesh.E
    B = mesh.B
    # NB: alias the current-density field as `J_field` rather than `j`,
    # because the loop variable `j` shadows it (and `B/E/J` are looked up
    # by name in the mesh struct, not by their field handles).
    J_field = mesh.j

    dx = config.Δ[1]
    dy = config.Δ[2]
    cdt_dx = dt / dx
    cdt_dy = dt / dy

    # Update E_x
    @inbounds for jj in 2:size(E, 2)-1, i in 2:size(E, 1)-1
        E[i, jj, 1] += cdt_dx * (B[i, jj+1, 3] - B[i, jj, 3])
        E[i, jj, 1] -= dt * J_field[i, jj, 1]
    end

    # Update E_y
    @inbounds for jj in 2:size(E, 2)-1, i in 2:size(E, 1)-1
        E[i, jj, 2] -= cdt_dy * (B[i+1, jj, 3] - B[i, jj, 3])
        E[i, jj, 2] -= dt * J_field[i, jj, 2]
    end

    # Update E_z
    @inbounds for jj in 2:size(E, 2)-1, i in 2:size(E, 1)-1
        E[i, jj, 3] += cdt_dx * (B[i+1, jj, 2] - B[i, jj, 2])
        E[i, jj, 3] -= cdt_dy * (B[i, jj+1, 1] - B[i, jj, 1])
        E[i, jj, 3] -= dt * J_field[i, jj, 3]
    end

    apply_boundary_conditions(mesh, :E)
end

function update_electric_field(mesh::AbstractMesh, fdtd::FDTD3D, dt::Float64)
    # 3D FDTD electric field update
    config = mesh.config
    E = mesh.E
    B = mesh.B
    J_field = mesh.j

    dx = config.Δ[1]
    dy = config.Δ[2]
    dz = config.Δ[3]
    cdt_dx = dt / dx
    cdt_dy = dt / dy
    cdt_dz = dt / dz

    # Update E_x
    @inbounds for k in 2:size(E, 3)-1, jj in 2:size(E, 2)-1, i in 2:size(E, 1)-1
        E[i, jj, k, 1] += cdt_dy * (B[i, jj+1, k, 3] - B[i, jj, k, 3])
        E[i, jj, k, 1] -= cdt_dz * (B[i, jj, k+1, 2] - B[i, jj, k, 2])
        E[i, jj, k, 1] -= dt * J_field[i, jj, k, 1]
    end

    # Update E_y
    @inbounds for k in 2:size(E, 3)-1, jj in 2:size(E, 2)-1, i in 2:size(E, 1)-1
        E[i, jj, k, 2] += cdt_dz * (B[i, jj, k+1, 1] - B[i, jj, k, 1])
        E[i, jj, k, 2] -= cdt_dx * (B[i+1, jj, k, 3] - B[i, jj, k, 3])
        E[i, jj, k, 2] -= dt * J_field[i, jj, k, 2]
    end

    # Update E_z
    @inbounds for k in 2:size(E, 3)-1, jj in 2:size(E, 2)-1, i in 2:size(E, 1)-1
        E[i, jj, k, 3] += cdt_dx * (B[i+1, jj, k, 2] - B[i, jj, k, 2])
        E[i, jj, k, 3] -= cdt_dy * (B[i, jj+1, k, 1] - B[i, jj, k, 1])
        E[i, jj, k, 3] -= dt * J_field[i, jj, k, 3]
    end

    apply_boundary_conditions(mesh, :E)
end

function update_magnetic_field(mesh::AbstractMesh, fdtd::FDTD1D, dt::Float64)
    # 1D FDTD magnetic field update
    config = mesh.config
    E = mesh.E
    B = mesh.B

    dx = config.Δ[1]
    cdt_dx = dt / dx

    # Update B_y / B_z components (1D)
    @inbounds for i in 2:size(B, 1)-1
        diff = E[i+1, 1] - E[i, 1]
        B[i, 2] -= cdt_dx * diff
        B[i, 3] += cdt_dx * diff
    end

    apply_boundary_conditions(mesh, :B)
end

function update_magnetic_field(mesh::AbstractMesh, fdtd::FDTD2D, dt::Float64)
    # 2D FDTD magnetic field update
    config = mesh.config
    E = mesh.E
    B = mesh.B

    dx = config.Δ[1]
    dy = config.Δ[2]
    cdt_dx = dt / dx
    cdt_dy = dt / dy

    # Update B_x
    @inbounds for j in 2:size(B, 2)-1, i in 2:size(B, 1)-1
        B[i, j, 1] -= cdt_dy * (E[i, j+1, 3] - E[i, j, 3])
    end

    # Update B_y
    @inbounds for j in 2:size(B, 2)-1, i in 2:size(B, 1)-1
        B[i, j, 2] += cdt_dx * (E[i+1, j, 3] - E[i, j, 3])
    end

    # Update B_z
    @inbounds for j in 2:size(B, 2)-1, i in 2:size(B, 1)-1
        B[i, j, 3] -= cdt_dx * (E[i+1, j, 2] - E[i, j, 2])
        B[i, j, 3] += cdt_dy * (E[i, j+1, 1] - E[i, j, 1])
    end

    apply_boundary_conditions(mesh, :B)
end

function update_magnetic_field(mesh::AbstractMesh, fdtd::FDTD3D, dt::Float64)
    # 3D FDTD magnetic field update
    config = mesh.config
    E = mesh.E
    B = mesh.B

    dx = config.Δ[1]
    dy = config.Δ[2]
    dz = config.Δ[3]
    cdt_dx = dt / dx
    cdt_dy = dt / dy
    cdt_dz = dt / dz

    # Update B_x
    @inbounds for k in 2:size(B, 3)-1, j in 2:size(B, 2)-1, i in 2:size(B, 1)-1
        B[i, j, k, 1] -= cdt_dy * (E[i, j+1, k, 3] - E[i, j, k, 3])
        B[i, j, k, 1] += cdt_dz * (E[i, j, k+1, 2] - E[i, j, k, 2])
    end

    # Update B_y
    @inbounds for k in 2:size(B, 3)-1, j in 2:size(B, 2)-1, i in 2:size(B, 1)-1
        B[i, j, k, 2] -= cdt_dz * (E[i, j, k+1, 1] - E[i, j, k, 1])
        B[i, j, k, 2] += cdt_dx * (E[i+1, j, k, 3] - E[i, j, k, 3])
    end

    # Update B_z
    @inbounds for k in 2:size(B, 3)-1, j in 2:size(B, 2)-1, i in 2:size(B, 1)-1
        B[i, j, k, 3] -= cdt_dx * (E[i+1, j, k, 2] - E[i, j, k, 2])
        B[i, j, k, 3] += cdt_dy * (E[i, j+1, k, 1] - E[i, j, k, 1])
    end

    apply_boundary_conditions(mesh, :B)
end

"""
$(TYPEDSIGNATURES)
Apply the configured boundary condition to a mesh field. The branch is
selected from the mesh configuration (`Periodic`, `Dirichlet`, `Vacuum`).
"""
function apply_boundary_conditions(mesh::AbstractMesh, field_name::Symbol)
    config = mesh.config
    boundary = config.boundary

    if boundary isa Periodic
        apply_periodic_boundary(mesh, field_name)
    elseif boundary isa Dirichlet
        apply_dirichlet_boundary(mesh, field_name)
    elseif boundary isa Vacuum
        apply_vacuum_boundary(mesh, field_name)
    end
end

"""
$(TYPEDSIGNATURES)
Apply periodic boundary conditions to a mesh field.
"""
function apply_periodic_boundary(mesh::AbstractMesh, field_name::Symbol)
    field = getproperty(mesh, field_name)
    field === nothing && return  # unallocated field — nothing to do

    nd = ndims(field.data)
    if nd == 1
        # 1D: data is array of length N with no extra vector dim
        @inbounds begin
            field.data[1] = field.data[end-1]
            field.data[end] = field.data[2]
        end
    elseif nd == 2
        # 1D vector field: data has shape (N, vec_dim)
        @inbounds begin
            field.data[1, :] .= field.data[end-1, :]
            field.data[end, :] .= field.data[2, :]
        end
    elseif nd == 3  # 2D field
        @inbounds begin
            # x boundaries
            field.data[1, :, :] .= field.data[end-1, :, :]
            field.data[end, :, :] .= field.data[2, :, :]

            # y boundaries
            field.data[:, 1, :] .= field.data[:, end-1, :]
            field.data[:, end, :] .= field.data[:, 2, :]
        end
    elseif nd == 4  # 3D field
        @inbounds begin
            # x boundaries
            field.data[1, :, :, :] .= field.data[end-1, :, :, :]
            field.data[end, :, :, :] .= field.data[2, :, :, :]

            # y boundaries
            field.data[:, 1, :, :] .= field.data[:, end-1, :, :]
            field.data[:, end, :, :] .= field.data[:, 2, :, :]

            # z boundaries
            field.data[:, :, 1, :] .= field.data[:, :, end-1, :]
            field.data[:, :, end, :] .= field.data[:, :, 2, :]
        end
    end
end

"""
$(TYPEDSIGNATURES)
Apply Dirichlet boundary conditions (zero field at the boundary) to a mesh
field.
"""
function apply_dirichlet_boundary(mesh::AbstractMesh, field_name::Symbol)
    field = getproperty(mesh, field_name)
    field === nothing && return  # unallocated field — nothing to do

    nd = ndims(field.data)
    z = zero(eltype(field.data))
    if nd == 1
        @inbounds begin
            field.data[1] = z
            field.data[end] = z
        end
    elseif nd == 2
        @inbounds begin
            field.data[1, :] .= z
            field.data[end, :] .= z
        end
    elseif nd == 3
        @inbounds begin
            field.data[1, :, :] .= z
            field.data[end, :, :] .= z
            field.data[:, 1, :] .= z
            field.data[:, end, :] .= z
        end
    elseif nd == 4
        @inbounds begin
            field.data[1, :, :, :] .= z
            field.data[end, :, :, :] .= z
            field.data[:, 1, :, :] .= z
            field.data[:, end, :, :] .= z
            field.data[:, :, 1, :] .= z
            field.data[:, :, end, :] .= z
        end
    end
end

"""
$(TYPEDSIGNATURES)
Apply vacuum (absorbing) boundary conditions to a mesh field by replicating
the nearest interior cells.
"""
function apply_vacuum_boundary(mesh::AbstractMesh, field_name::Symbol)
    field = getproperty(mesh, field_name)
    field === nothing && return  # unallocated field — nothing to do

    nd = ndims(field.data)
    if nd == 1
        @inbounds begin
            field.data[1] = field.data[2]
            field.data[end] = field.data[end-1]
        end
    elseif nd == 2
        @inbounds begin
            field.data[1, :] .= field.data[2, :]
            field.data[end, :] .= field.data[end-1, :]
        end
    elseif nd == 3
        @inbounds begin
            field.data[1, :, :] .= field.data[2, :, :]
            field.data[end, :, :] .= field.data[end-1, :, :]
            field.data[:, 1, :] .= field.data[:, 2, :]
            field.data[:, end, :] .= field.data[:, end-1, :]
        end
    elseif nd == 4
        @inbounds begin
            field.data[1, :, :, :] .= field.data[2, :, :, :]
            field.data[end, :, :, :] .= field.data[end-1, :, :, :]
            field.data[:, 1, :, :] .= field.data[:, 2, :, :]
            field.data[:, end, :, :] .= field.data[:, end-1, :, :]
            field.data[:, :, 1, :] .= field.data[:, :, 2, :]
            field.data[:, :, end, :] .= field.data[:, :, end-1, :]
        end
    end
end

"""
$(TYPEDSIGNATURES)
Run an FDTD simulation for the given number of steps with fixed timestep.
"""
function fdtd_simulate(mesh::AbstractMesh, fdtd::AbstractFDTD, dt::Float64, steps::Int)
    for i in 1:steps
        fdtd_step(mesh, fdtd, dt)
    end
end
