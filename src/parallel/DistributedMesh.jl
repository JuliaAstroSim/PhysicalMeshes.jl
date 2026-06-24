struct DistributedMesh{I, VI, V, U, D, M <: AbstractMesh{U}} <: AbstractMesh{U}
    config::MeshConfig{I,VI,V,U}
    local_mesh::M
    rank::Int
    size::Int
    neighbors::Vector{Int}

    # Explicit inner constructor.  The `D` type parameter (device kind
    # of the local config) is not directly bound by any field type — the
    # field `config` uses `MeshConfig{I,VI,V,U}` which leaves `D` as a
    # partial-instantiation parameter — so Julia's auto-generated
    # default positional constructor cannot always infer `D`, and calls
    # like `DistributedMesh(cfg, lm, r, s, n)` raise MethodError.
    # Pattern-matching the full `MeshConfig{I,VI,V,U,D}` here lets us
    # bind `D` from the actual `config` argument and reconstruct the
    # struct unambiguously.
    function DistributedMesh(
        config::MeshConfig{I,VI,V,U,D},
        local_mesh::M,
        rank::Int,
        size::Int,
        neighbors::Vector{Int},
    ) where {I,VI,V,U,D,M<:AbstractMesh{U}}
        return new{I,VI,V,U,D,M}(config, local_mesh, rank, size, neighbors)
    end
end

"""
$(TYPEDSIGNATURES)
Construct a `DistributedMesh` whose local subdomain corresponds to
the calling rank. The constructor relies on `calculate_local_config`
and `partition_neighbors` to derive the local extents and the rank's
neighbours from the global config.

`kw...` is forwarded to the underlying `MeshCartesianStatic` (e.g.,
`mode`, `assignment`, `boundary`).
"""
function DistributedMesh(config::MeshConfig, units = nothing; kw...)
    # 0-based rank, size includes the master process.
    rank = myid() - 1
    size = nworkers() + 1

    local_config = calculate_local_config(config, rank, size)
    local_mesh   = MeshCartesianStatic(local_config, nothing, units; kw...)
    neighbors    = partition_neighbors(rank, size, config.dim)

    return DistributedMesh(config, local_mesh, rank, size, neighbors)
end

function calculate_local_config(config::MeshConfig, rank::Int, size::Int)
    # Simple 1D decomposition for now
    dim = config.dim
    N = config.N
    Min = config.Min
    Max = config.Max
    Δ = config.Δ

    # Calculate local size and offset
    if dim == 1
        local_Nx = N[1] ÷ size
        remainder = N[1] % size
        local_Nx += rank < remainder ? 1 : 0
        offset_x = rank < remainder ? rank * (local_Nx) : rank * (local_Nx) + remainder

        local_Min = SVector(Min[1] + offset_x * Δ[1])
        local_Max = SVector(local_Min[1] + local_Nx * Δ[1])
        local_N = SVector(local_Nx)
    elseif dim == 2
        # Simple row-wise decomposition
        rows = Int(ceil(sqrt(size)))
        cols = Int(ceil(size / rows))
        row = rank ÷ cols
        col = rank % cols

        local_Nx = N[1] ÷ cols
        local_Ny = N[2] ÷ rows
        remainder_x = N[1] % cols
        remainder_y = N[2] % rows

        local_Nx += col < remainder_x ? 1 : 0
        local_Ny += row < remainder_y ? 1 : 0

        offset_x = col < remainder_x ? col * (local_Nx) : col * (local_Nx) + remainder_x
        offset_y = row < remainder_y ? row * (local_Ny) : row * (local_Ny) + remainder_y

        local_Min = SVector(Min[1] + offset_x * Δ[1], Min[2] + offset_y * Δ[2])
        local_Max = SVector(local_Min[1] + local_Nx * Δ[1], local_Min[2] + local_Ny * Δ[2])
        local_N = SVector(local_Nx, local_Ny)
    else # dim == 3
        # Simple 3D decomposition
        # This is a basic implementation, can be optimized
        local_Nx = N[1] ÷ size
        local_Ny = N[2]
        local_Nz = N[3]
        remainder = N[1] % size
        local_Nx += rank < remainder ? 1 : 0
        offset_x = rank < remainder ? rank * (local_Nx) : rank * (local_Nx) + remainder

        local_Min = SVector(Min[1] + offset_x * Δ[1], Min[2], Min[3])
        local_Max = SVector(local_Min[1] + local_Nx * Δ[1], Max[2], Max[3])
        local_N = SVector(local_Nx, local_Ny, local_Nz)
    end

    # Create local config
    return MeshConfig(
        config.mode,
        config.assignment,
        config.boundary,
        config.enlarge,
        config.device,
        config.units,
        config.dim,
        config.NG,
        Δ,
        local_Min,
        local_Max,
        local_N,
        local_N .+ 2 * config.NG
    )
end

"""
$(TYPEDSIGNATURES)
Return the rank neighbours of `rank` in a partition of size `size` for
a `dim`-dimensional grid. 1D and 3D partitions use a linear row layout;
2D partitions use the same `ceil(sqrt(size))` row / `ceil(size/rows)`
column layout adopted by `calculate_local_config`.
"""
function partition_neighbors(rank::Int, size::Int, dim::Int)
    neighbors = Int[]

    if dim == 1
        # Left neighbor
        if rank > 0
            push!(neighbors, rank - 1)
        end
        # Right neighbor
        if rank < size - 1
            push!(neighbors, rank + 1)
        end
    elseif dim == 2
        rows = Int(ceil(sqrt(size)))
        cols = Int(ceil(size / rows))
        row = rank ÷ cols
        col = rank % cols

        # Top neighbor
        if row > 0
            neighbor = (row - 1) * cols + col
            if 0 <= neighbor < size
                push!(neighbors, neighbor)
            end
        end
        # Bottom neighbor
        if row < rows - 1
            neighbor = (row + 1) * cols + col
            if 0 <= neighbor < size
                push!(neighbors, neighbor)
            end
        end
        # Left neighbor
        if col > 0
            neighbor = row * cols + (col - 1)
            0 <= neighbor < size && push!(neighbors, neighbor)
        end
        # Right neighbor
        if col < cols - 1
            neighbor = row * cols + (col + 1)
            0 <= neighbor < size && push!(neighbors, neighbor)
        end
    else # dim == 3
        # Simple 1D decomposition neighbors
        if rank > 0
            push!(neighbors, rank - 1)
        end
        if rank < size - 1
            push!(neighbors, rank + 1)
        end
    end

    return filter(n -> 0 <= n < size, neighbors)
end

function Base.show(io::IO, mesh::DistributedMesh)
    print(io, "DistributedMesh on rank $(mesh.rank)/$(mesh.size)\n")
    print(io, "Local mesh: $(mesh.local_mesh)\n")
    print(io, "Neighbors: $(mesh.neighbors)")
end

# Field accessors: use getproperty so that `mesh.field` delegates to the
# underlying local mesh for the common FDTD/PIC use cases.
@inline Base.getproperty(mesh::DistributedMesh, name::Symbol) =
    name in fieldnames(typeof(mesh)) ? Base.getfield(mesh, name) :
    Base.getproperty(mesh.local_mesh, name)

# AbstractMesh fallbacks
number_of_nodes(mesh::DistributedMesh) = number_of_nodes(mesh.local_mesh)
number_of_cells(mesh::DistributedMesh) = number_of_cells(mesh.local_mesh)
number_of_fields(mesh::DistributedMesh) = number_of_fields(mesh.local_mesh)
