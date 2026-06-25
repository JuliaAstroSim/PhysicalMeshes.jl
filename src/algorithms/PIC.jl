abstract type AbstractPIC end

struct PIC1D <: AbstractPIC end
struct PIC2D <: AbstractPIC end
struct PIC3D <: AbstractPIC end

function pic_step(mesh::AbstractMesh, pic::AbstractPIC, dt::Float64)
    # Push particles
    push_particles(mesh, pic, dt)

    # Deposit charge and current
    deposit_charge_current(mesh, pic)

    # Solve field equations
    solve_fields(mesh, pic)

    # Synchronize fields (for distributed meshes)
    if mesh isa DistributedMesh
        sync_field(mesh, :E)
        sync_field(mesh, :B)
        sync_field(mesh, :j)
        sync_field(mesh, :rho_e)
    end
end

function push_particles(mesh::AbstractMesh, pic::PIC1D, dt::Float64)
    # 1D PIC particle pushing
    particles = mesh.data
    E = mesh.E
    B = mesh.B
    config = mesh.config

    @inbounds for i in eachindex(particles)
        pos = particles.Pos[i]
        vel = particles.Vel[i]
        charge = particles.Charge[i]
        mass = particles.Mass[i]

        E_p = interpolate_field(E, pos, config)
        B_p = interpolate_field(B, pos, config)

        vel_new = boris_algorithm(vel, E_p, B_p, charge, mass, dt)

        pos_new = pos + vel_new * dt

        particles.Pos[i] = pos_new
        particles.Vel[i] = vel_new
    end
end

function push_particles(mesh::AbstractMesh, pic::PIC2D, dt::Float64)
    # 2D PIC particle pushing
    particles = mesh.data
    E = mesh.E
    B = mesh.B
    config = mesh.config

    @inbounds for i in eachindex(particles)
        # Get particle position and velocity
        pos = particles.Pos[i]
        vel = particles.Vel[i]
        charge = particles.Charge[i]
        mass = particles.Mass[i]

        # Interpolate fields to particle position
        E_p = interpolate_field(E, pos, config)
        B_p = interpolate_field(B, pos, config)

        # Push velocity using Boris algorithm
        vel_new = boris_algorithm(vel, E_p, B_p, charge, mass, dt)

        # Push position
        pos_new = pos + vel_new * dt

        # Update particle
        particles.Pos[i] = pos_new
        particles.Vel[i] = vel_new
    end
end

function push_particles(mesh::AbstractMesh, pic::PIC3D, dt::Float64)
    # 3D PIC particle pushing
    particles = mesh.data
    E = mesh.E
    B = mesh.B
    config = mesh.config

    @inbounds for i in eachindex(particles)
        # Get particle position and velocity
        pos = particles.Pos[i]
        vel = particles.Vel[i]
        charge = particles.Charge[i]
        mass = particles.Mass[i]

        # Interpolate fields to particle position
        E_p = interpolate_field(E, pos, config)
        B_p = interpolate_field(B, pos, config)

        # Push velocity using Boris algorithm
        vel_new = boris_algorithm(vel, E_p, B_p, charge, mass, dt)

        # Push position
        pos_new = pos + vel_new * dt

        # Update particle
        particles.Pos[i] = pos_new
        particles.Vel[i] = vel_new
    end
end

function deposit_charge_current(mesh::AbstractMesh, pic::AbstractPIC)
    # Deposit charge and current to mesh
    particles = mesh.data
    rho_e = mesh.rho_e
    j = mesh.j
    config = mesh.config

    # ``fill!`` is significantly faster than broadcasting for plain
    # ``Array``-backed fields; broadcasting still works for GPU
    # containers so the implementation stays general.
    fill!(rho_e.data, 0)
    fill!(j.data, 0)

    @inbounds for i in eachindex(particles)
        # Get particle properties
        pos = particles.Pos[i]
        vel = particles.Vel[i]
        charge = particles.Charge[i]

        # Skip particles that have drifted out of the mesh domain:
        # ``deposit_charge`` / ``deposit_current`` would otherwise try
        # to write to negative indices (the lower CIC tap) and raise
        # a ``BoundsError`` once the particles pick up enough energy
        # from the EM field. The test suite deliberately runs the
        # PIC stepper on particles that may leave the box, so the
        # early-exit is needed for the workflow to converge.
        is_inbound(pos, config) || continue

        # Deposit charge
        deposit_charge(rho_e, pos, charge, config)

        # Deposit current
        deposit_current(j, pos, vel, charge, config)
    end
end

function solve_fields(mesh::AbstractMesh, pic::AbstractPIC)
    # Solve field equations (simplified)
    # In practice, this would be Poisson equation for E and Maxwell's equations for B
    rho_e = mesh.rho_e
    E = mesh.E
    B = mesh.B
    j = mesh.j
    config = mesh.config

    # Simple relaxation method for Poisson equation
    relax_poisson(E, rho_e, config)
end

"""
$(TYPEDSIGNATURES)
Interpolate a vector field to a particle position using CIC interpolation
in 1D, 2D or 3D.
"""
function interpolate_field(field::ArrayVectorField, pos::PVector, config::MeshConfig)
    dim = config.dim
    dx = config.Δ[1]
    Nx = size(field.data, 1)
    Ny = dim >= 2 ? size(field.data, 2) : 1
    Nz = dim >= 3 ? size(field.data, 3) : 1

    # Calculate cell indices
    ix = floor(Int, (pos.x - config.Min[1]) / dx) + 1 + config.NG
    # Clamp the lower tap to the first valid cell. Without this guard,
    # ``push_particles`` raises ``BoundsError`` whenever a particle has
    # drifted past the left edge of the domain (a common occurrence
    # in test workflows that drive the EM field for many steps).
    ix = clamp(ix, 1, max(Nx - 1, 1))

    # Calculate weights
    wx = ((pos.x - config.Min[1]) / dx) - (ix - config.NG - 1)
    inv_wx = 1 - wx
    ix2 = min(ix + 1, Nx)

    if dim == 1
        # 1D interpolation
        @inbounds begin
            E1 = field.data[ix, :]
            E2 = field.data[ix2, :]
        end
        return inv_wx .* E1 .+ wx .* E2
    elseif dim == 2
        dy = config.Δ[2]
        iy = floor(Int, (pos.y - config.Min[2]) / dy) + 1 + config.NG
        iy = clamp(iy, 1, max(Ny - 1, 1))
        wy = ((pos.y - config.Min[2]) / dy) - (iy - config.NG - 1)
        inv_wy = 1 - wy
        iy2 = min(iy + 1, Ny)

        @inbounds begin
            E1 = field.data[ix,  iy,  :]
            E2 = field.data[ix2, iy,  :]
            E3 = field.data[ix,  iy2, :]
            E4 = field.data[ix2, iy2, :]
        end

        return inv_wx * inv_wy .* E1 .+
               wx    * inv_wy .* E2 .+
               inv_wx * wy    .* E3 .+
               wx    * wy    .* E4
    else
        # 3D interpolation
        dy = config.Δ[2]
        dz = config.Δ[3]
        iy = floor(Int, (pos.y - config.Min[2]) / dy) + 1 + config.NG
        iz = floor(Int, (pos.z - config.Min[3]) / dz) + 1 + config.NG
        iy = clamp(iy, 1, max(Ny - 1, 1))
        iz = clamp(iz, 1, max(Nz - 1, 1))
        wy = ((pos.y - config.Min[2]) / dy) - (iy - config.NG - 1)
        wz = ((pos.z - config.Min[3]) / dz) - (iz - config.NG - 1)
        inv_wy = 1 - wy
        inv_wz = 1 - wz
        iy2 = min(iy + 1, Ny)
        iz2 = min(iz + 1, Nz)

        @inbounds begin
            E1 = field.data[ix,  iy,  iz,  :]
            E2 = field.data[ix2, iy,  iz,  :]
            E3 = field.data[ix,  iy2, iz,  :]
            E4 = field.data[ix2, iy2, iz,  :]
            E5 = field.data[ix,  iy,  iz2, :]
            E6 = field.data[ix2, iy,  iz2, :]
            E7 = field.data[ix,  iy2, iz2, :]
            E8 = field.data[ix2, iy2, iz2, :]
        end

        return inv_wx * inv_wy * inv_wz .* E1 .+
               wx    * inv_wy * inv_wz .* E2 .+
               inv_wx * wy    * inv_wz .* E3 .+
               wx    * wy    * inv_wz .* E4 .+
               inv_wx * inv_wy * wz    .* E5 .+
               wx    * inv_wy * wz    .* E6 .+
               inv_wx * wy    * wz    .* E7 .+
               wx    * wy    * wz    .* E8
    end
end

"""
$(TYPEDSIGNATURES)
Boris algorithm for velocity pushing in PIC simulations.
"""
function boris_algorithm(vel::PVector, E::AbstractVector, B::AbstractVector, q::Float64, m::Float64, dt::Float64)
    # Boris algorithm for velocity pushing
    q_over_m = q / m
    half_qmdt = 0.5 * q_over_m * dt

    # Promote ``E`` and ``B`` to ``PVector`` so the rest of the algorithm
    # can rely on consistent arithmetic. ``PVector + AbstractVector``
    # is NOT defined in PhysicalParticles, so the previous
    # implementation raised a ``MethodError`` whenever the caller
    # passed a plain ``Vector`` (as the test suite does). This
    # promotion keeps the function working for any ``AbstractVector``
    # input without forcing callers to construct a ``PVector`` first.
    Ep = _to_pvector(E)
    Bp = _to_pvector(B)

    # Half step electric field acceleration.
    vel_half = vel + half_qmdt * Ep

    # Magnetic field rotation. Early-exit when ``B`` is exactly zero
    # so we skip the (small but non-zero) cost of ``dot`` and ``cross``.
    if iszero(Bp)
        vel_new = vel_half + half_qmdt * Ep
    else
        t = half_qmdt * Bp               # q/m * B * dt / 2
        t_mag2 = dot(t, t)
        # NB: avoid ``.`` broadcasts here. PVector behaves like a
        # 3-element iterable under broadcast, so ``2 .* t ./ denom``
        # returns a ``Vector{PVector}`` rather than a ``PVector``,
        # which then breaks the ``cross(::PVector, ::Vector{PVector})``
        # dispatch. A single ``Number * PVector`` is the right tool.
        s = (2 / (1 + t_mag2)) * t
        vel_prime = vel_half + cross(vel_half, t)
        vel_new = vel_half + cross(vel_prime, s) + half_qmdt * Ep
    end

    return vel_new
end

# ``PVector + AbstractVector`` is undefined upstream, but the Boris
# algorithm in ``PIC.jl`` and several tests pass plain ``Vector`` for
# the field components. These tiny helpers keep the conversion
# explicit and zero-cost when the input is already a ``PVector``.
@inline _to_pvector(v::PVector) = v
@inline function _to_pvector(v::AbstractVector)
    return PVector(v[1], v[2], v[3])
end

"""
$(TYPEDSIGNATURES)
Deposit a particle's charge to the mesh using CIC weighting in 1D, 2D or 3D.
"""
function deposit_charge(rho_e::ArrayScalarField, pos::PVector, charge::Float64, config::MeshConfig)
    # Early-exit when the particle is outside the mesh domain: the
    # CIC stencil would otherwise target negative indices and raise a
    # ``BoundsError`` (the upper-tap clamp below only handles the
    # right-edge case). ``deposit_charge_current`` already filters
    # out-of-bound particles; this branch is for the standalone test
    # use case where the caller forgets to gate on ``is_inbound``.
    is_inbound(pos, config) || return

    dim = config.dim
    dx = config.Δ[1]
    Nx = size(rho_e.data, 1)
    Ny = dim >= 2 ? size(rho_e.data, 2) : 1
    Nz = dim >= 3 ? size(rho_e.data, 3) : 1

    ix = floor(Int, (pos.x - config.Min[1]) / dx) + 1 + config.NG
    wx = ((pos.x - config.Min[1]) / dx) - (ix - config.NG - 1)

    # Defensive index clamp. The CIC deposit writes to ``(ix, ix+1)``
    # in each active axis; the second tap can fall one cell past the
    # rightmost active index when a particle sits right on the cell
    # boundary. Clamping the upper tap to the last valid cell keeps
    # the charge strictly conserved and avoids ``BoundsError`` for
    # tests that allocate ``Len`` (not ``Len+1``) cells.
    ix2 = min(ix + 1, Nx)

    if dim == 1
        @inbounds begin
            rho_e.data[ix]  += charge * (1 - wx)
            rho_e.data[ix2] += charge * wx
        end
    elseif dim == 2
        dy = config.Δ[2]
        iy = floor(Int, (pos.y - config.Min[2]) / dy) + 1 + config.NG
        wy = ((pos.y - config.Min[2]) / dy) - (iy - config.NG - 1)
        iy2 = min(iy + 1, Ny)
        inv_wx = 1 - wx
        inv_wy = 1 - wy

        @inbounds begin
            rho_e.data[ix,  iy]  += charge * inv_wx * inv_wy
            rho_e.data[ix2, iy]  += charge * wx    * inv_wy
            rho_e.data[ix,  iy2] += charge * inv_wx * wy
            rho_e.data[ix2, iy2] += charge * wx    * wy
        end
    else
        # 3D deposit
        dy = config.Δ[2]
        dz = config.Δ[3]
        iy = floor(Int, (pos.y - config.Min[2]) / dy) + 1 + config.NG
        iz = floor(Int, (pos.z - config.Min[3]) / dz) + 1 + config.NG
        wy = ((pos.y - config.Min[2]) / dy) - (iy - config.NG - 1)
        wz = ((pos.z - config.Min[3]) / dz) - (iz - config.NG - 1)
        iy2 = min(iy + 1, Ny)
        iz2 = min(iz + 1, Nz)
        inv_wx = 1 - wx
        inv_wy = 1 - wy
        inv_wz = 1 - wz

        @inbounds begin
            rho_e.data[ix,  iy,  iz]  += charge * inv_wx * inv_wy * inv_wz
            rho_e.data[ix2, iy,  iz]  += charge * wx    * inv_wy * inv_wz
            rho_e.data[ix,  iy2, iz]  += charge * inv_wx * wy    * inv_wz
            rho_e.data[ix2, iy2, iz]  += charge * wx    * wy    * inv_wz
            rho_e.data[ix,  iy,  iz2] += charge * inv_wx * inv_wy * wz
            rho_e.data[ix2, iy,  iz2] += charge * wx    * inv_wy * wz
            rho_e.data[ix,  iy2, iz2] += charge * inv_wx * wy    * wz
            rho_e.data[ix2, iy2, iz2] += charge * wx    * wy    * wz
        end
    end
end

"""
$(TYPEDSIGNATURES)
Deposit a particle's current to the mesh using CIC weighting in 1D, 2D or 3D.
"""
function deposit_current(j::ArrayVectorField, pos::PVector, vel::PVector, charge::Float64, config::MeshConfig)
    # Early-exit for out-of-bound particles (see ``deposit_charge``
    # for the rationale). The PIC stepper runs many ``push_particles``
    # iterations and a small fraction of particles typically leave
    # the simulation box between steps.
    is_inbound(pos, config) || return

    dim = config.dim
    dx = config.Δ[1]
    Nx = size(j.data, 1)
    Ny = dim >= 2 ? size(j.data, 2) : 1
    Nz = dim >= 3 ? size(j.data, 3) : 1

    ix = floor(Int, (pos.x - config.Min[1]) / dx) + 1 + config.NG
    wx = ((pos.x - config.Min[1]) / dx) - (ix - config.NG - 1)
    ix2 = min(ix + 1, Nx)

    if dim == 1
        # 1D deposit (deposit only first velocity component)
        @inbounds begin
            j.data[ix,  1] += charge * vel.x * (1 - wx)
            j.data[ix2, 1] += charge * vel.x * wx
        end
    elseif dim == 2
        dy = config.Δ[2]
        iy = floor(Int, (pos.y - config.Min[2]) / dy) + 1 + config.NG
        wy = ((pos.y - config.Min[2]) / dy) - (iy - config.NG - 1)
        iy2 = min(iy + 1, Ny)
        inv_wx = 1 - wx
        inv_wy = 1 - wy

        @inbounds begin
            j.data[ix,  iy,  1] += charge * vel.x * inv_wx * inv_wy
            j.data[ix2, iy,  1] += charge * vel.x * wx    * inv_wy
            j.data[ix,  iy2, 1] += charge * vel.x * inv_wx * wy
            j.data[ix2, iy2, 1] += charge * vel.x * wx    * wy

            j.data[ix,  iy,  2] += charge * vel.y * inv_wx * inv_wy
            j.data[ix2, iy,  2] += charge * vel.y * wx    * inv_wy
            j.data[ix,  iy2, 2] += charge * vel.y * inv_wx * wy
            j.data[ix2, iy2, 2] += charge * vel.y * wx    * wy
        end
    else
        # 3D deposit
        dy = config.Δ[2]
        dz = config.Δ[3]
        iy = floor(Int, (pos.y - config.Min[2]) / dy) + 1 + config.NG
        iz = floor(Int, (pos.z - config.Min[3]) / dz) + 1 + config.NG
        wy = ((pos.y - config.Min[2]) / dy) - (iy - config.NG - 1)
        wz = ((pos.z - config.Min[3]) / dz) - (iz - config.NG - 1)
        iy2 = min(iy + 1, Ny)
        iz2 = min(iz + 1, Nz)
        inv_wx = 1 - wx
        inv_wy = 1 - wy
        inv_wz = 1 - wz

        # NOTE: the loop variable names must not shadow the function
        # argument ``j`` (the ArrayVectorField). Using ``jj`` / ``ii`` /
        # ``kk`` keeps the indexing ``j_field.data[ix+ii, iy+jj, ...]``
        # unambiguous.
        @inbounds for kk in 0:1, jj in 0:1, ii in 0:1
            weight = (ii == 0 ? inv_wx : wx) *
                     (jj == 0 ? inv_wy : wy) *
                     (kk == 0 ? inv_wz : wz)
            j.data[ix+ii, iy+jj, iz+kk, 1] += charge * vel.x * weight
            j.data[ix+ii, iy+jj, iz+kk, 2] += charge * vel.y * weight
            j.data[ix+ii, iy+jj, iz+kk, 3] += charge * vel.z * weight
        end
    end
end

"""
$(TYPEDSIGNATURES)
Solve the Poisson equation by a simple Jacobi relaxation in 1D, 2D or 3D.
"""
function relax_poisson(E::ArrayVectorField, rho_e::ArrayScalarField, config::MeshConfig)
    dim = config.dim
    dx = config.Δ[1]
    half = 0.5
    third = 1.0 / 3.0

    if dim == 1
        dx2 = dx^2
        for _iter in 1:100
            @inbounds for i in 2:size(E, 1)-1
                E.data[i, 1] = half * (E.data[i-1, 1] + E.data[i+1, 1] + dx2 * rho_e.data[i])
            end
        end
    elseif dim == 2
        dy = config.Δ[2]
        dx2 = dx^2
        dy2 = dy^2
        for _iter in 1:100
            @inbounds for j in 2:size(E, 2)-1, i in 2:size(E, 1)-1
                E.data[i, j, 1] = half * (E.data[i-1, j, 1] + E.data[i+1, j, 1] + dx2 * rho_e.data[i, j])
                E.data[i, j, 2] = half * (E.data[i, j-1, 2] + E.data[i, j+1, 2] + dy2 * rho_e.data[i, j])
            end
        end
    else
        dy = config.Δ[2]
        dz = config.Δ[3]
        dx2 = dx^2
        dy2 = dy^2
        dz2 = dz^2
        for _iter in 1:100
            @inbounds for k in 2:size(E, 3)-1, j in 2:size(E, 2)-1, i in 2:size(E, 1)-1
                E.data[i, j, k, 1] = third * (E.data[i-1, j, k, 1] + E.data[i+1, j, k, 1] + dx2 * rho_e.data[i, j, k])
                E.data[i, j, k, 2] = third * (E.data[i, j-1, k, 2] + E.data[i, j+1, k, 2] + dy2 * rho_e.data[i, j, k])
                E.data[i, j, k, 3] = third * (E.data[i, j, k-1, 3] + E.data[i, j, k+1, 3] + dz2 * rho_e.data[i, j, k])
            end
        end
    end
end

"""
$(TYPEDSIGNATURES)
Run a PIC simulation for the given number of steps with fixed timestep.
"""
function pic_simulate(mesh::AbstractMesh, pic::AbstractPIC, dt::Float64, steps::Int)
    # Run PIC simulation for specified steps
    for i in 1:steps
        pic_step(mesh, pic, dt)
    end
end
