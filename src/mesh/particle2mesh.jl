function NGPinfo(meshpos, config, pos::AbstractVector{T}) where T<:Number
    id = round.(Int, (pos .- config.Min) ./ config.Δ) .+ 1 .+ config.NG
    return id
end

function CICinfo(meshpos, config, pos::AbstractVector{T}) where T<:Number
    # Here, left means smaller coordinate value in that dimension
    # find the index of nearest left and right vertex
    idl = floor.(Int, (pos .- config.Min) ./ config.Δ) .+ 1 .+ config.NG
    idr = idl .+ 1

    pl = SVector(meshpos[idl...])       # position of left vertex

    rl = (pl .- pos) ./ config.Δ .+ 1.0  # assignment ratio of left vertex
    rr = 1.0 .- rl                       # assignment ratio of right vertex

    id = [idl, idr]                      # pack info into array
    r = [rl, rr]

    return id, r
end

function TSCinfo(meshpos, config, pos::AbstractVector{T}) where T<:Number
    idl = floor.(Int, (pos .- config.Min) ./ config.Δ .- 0.5) .+ 1 .+ config.NG  #TODO check -0.5
    idm = idl .+ 1
    idr = idm .+ 1

    pl = SVector(meshpos[idl...])

    rl = (pl .- pos .+ 1.5 .* config.Δ) .^ 2 ./ 2 ./ config.Δ ./ config.Δ
    rr = (pl .- pos .+ 0.5 .* config.Δ) .^ 2 ./ 2 ./ config.Δ ./ config.Δ
    rm = 1.0 .- rl .- rr    # assignment ratio of middle vertex

    id = [idl, idm, idr]
    r = [rl, rm, rr]
    return id, r
end

"""
$(TYPEDSIGNATURES)
Determine whether a 3D position is inside the mesh domain.
"""
@inline function is_inbound(pos::PVector, config::MeshConfig)
    # Iterate over the mesh's actual dimensionality. The previous version
    # unconditionally read `config.Max[1..3]`, which BoundsError'd for a
    # 1D or 2D `MeshConfig` (where `Max` has length 1 or 2).
    @inbounds for i in eachindex(config.Min)
        coord = i == 1 ? pos.x : i == 2 ? pos.y : pos.z
        if coord > config.Max[i] || coord < config.Min[i]
            return false
        end
    end
    return true
end

function outbound_list(pos::AbstractArray, m::MeshCartesianStatic)
    config = m.config
    list = Int[]
    for i in eachindex(pos)
        if !is_inbound(pos[i], config)
            push!(list, i)
        end
    end
    return list
end

outbound_list(m::MeshCartesianStatic) = CUDA.@allowscalar outbound_list(Array(m.data.Pos), m)

@inline _mesh_index(id, k) = id[k]

function particle2mesh!(meshpos, config, pos::AbstractVector{T}, meshdata, particledata, ::VertexMode, ::NGP) where T<:Number
    id = NGPinfo(meshpos, config, pos)
    # Use only the first `config.dim` entries; data may be a higher-dim array
    @inbounds meshdata[ntuple(k -> _mesh_index(id, k), config.dim)...] += particledata
end

function particle2mesh!(meshpos, config, pos::AbstractVector{T}, meshdata, particledata, ::VertexMode, ::CIC) where T<:Number
    id, r = CICinfo(meshpos, config, pos)
    _scatter_cic(meshdata, id, r, particledata, Val(config.dim))
end

function particle2mesh!(meshpos, config, pos::AbstractVector{T}, meshdata, particledata, ::VertexMode, ::TSC) where T<:Number
    id, r = TSCinfo(meshpos, config, pos)
    _scatter_tsc(meshdata, id, r, particledata, Val(config.dim))
end

# Dimension-generic scatter helpers
@inline function _scatter_cic(meshdata, id, r, particledata, ::Val{1})
    @inbounds for i in 1:2
        meshdata[id[i][1]] += particledata * r[i][1]
    end
end

@inline function _scatter_cic(meshdata, id, r, particledata, ::Val{2})
    @inbounds for j in 1:2, i in 1:2
        meshdata[id[i][1], id[j][2]] += particledata * r[i][1] * r[j][2]
    end
end

@inline function _scatter_cic(meshdata, id, r, particledata, ::Val{3})
    @inbounds for k in 1:2, j in 1:2, i in 1:2
        meshdata[id[i][1], id[j][2], id[k][3]] += particledata * r[i][1] * r[j][2] * r[k][3]
    end
end

@inline function _scatter_tsc(meshdata, id, r, particledata, ::Val{1})
    @inbounds for i in 1:3
        meshdata[id[i][1]] += particledata * r[i][1]
    end
end

@inline function _scatter_tsc(meshdata, id, r, particledata, ::Val{2})
    @inbounds for j in 1:3, i in 1:3
        meshdata[id[i][1], id[j][2]] += particledata * r[i][1] * r[j][2]
    end
end

@inline function _scatter_tsc(meshdata, id, r, particledata, ::Val{3})
    @inbounds for k in 1:3, j in 1:3, i in 1:3
        meshdata[id[i][1], id[j][2], id[k][3]] += particledata * r[i][1] * r[j][2] * r[k][3]
    end
end

"""
$(TYPEDSIGNATURES)

Assign mesh data by mass assignment shemes

## Example
```julia
assignmesh(data, m, :Mass, :rho)
```
"""
Base.@propagate_inbounds function assignmesh(particles::StructArray, mesh::MeshCartesianStatic, symbolParticle::Symbol, symbolMesh::Symbol)
    config = mesh.config
    meshfield = getproperty(mesh, symbolMesh)

    # Zero out the field, supporting any AbstractField subtype
    _zero_mesh_field!(meshfield)

    CUDA.@allowscalar meshpos = Array(mesh.pos)
    CUDA.@allowscalar meshdata = Array(meshfield.data)
    CUDA.@allowscalar particlepos = Array(particles.Pos)
    CUDA.@allowscalar particledata = Array(getproperty(particles, symbolParticle))

    # Pre-compute the cell-volume normalisation factor once, outside the
    # particle loop. ``prod(config.Δ)`` was previously re-evaluated per
    # particle, which the JIT can normally hoist — but only when the
    # shape and type are known statically. Hoisting it explicitly here
    # guarantees it and documents the intent.
    inv_dV = one(eltype(particledata)) / prod(config.Δ)

    @inbounds for i in eachindex(particles)
        rho = particledata[i] * inv_dV
        if is_inbound(particlepos[i], config)
            pos = SVector(particlepos[i])
            particle2mesh!(meshpos, config, pos, meshdata, rho, config.mode, config.assignment)
        end
    end

    if meshfield.data isa CuArray
        meshfield.data .= cu(meshdata)
    else
        meshfield.data .= meshdata
    end
end

assignmesh(m::MeshCartesianStatic) = assignmesh(m.data, m, :Mass, :rho)

"""
    _zero_mesh_field!(field)

Reset every supported mesh field to zero. Dispatches on the concrete
`AbstractField` subtype so both scalar and vector fields are handled.
"""
@inline _zero_mesh_field!(field::ArrayScalarField) = (field.data .= zero(eltype(field.data)); nothing)
@inline _zero_mesh_field!(field::ArrayVectorField) = (field.data .= zero(eltype(field.data)); nothing)
@inline _zero_mesh_field!(field::ArrayTensorField) = (field.data .= zero(eltype(field.data)); nothing)
@inline _zero_mesh_field!(::Nothing) = nothing

function mesh2particle(meshpos, config, meshdata, pos::AbstractVector{T}, ::VertexMode, ::NGP) where T<:Number
    id = NGPinfo(meshpos, config, pos)
    return meshdata[ntuple(k -> id[k], config.dim)...]
end

function mesh2particle(meshpos, config, meshdata, pos::AbstractVector{T}, ::VertexMode, ::CIC) where T<:Number
    id, r = CICinfo(meshpos, config, pos)
    return _gather(meshdata, id, r, Val(config.dim), Val(2))
end

function mesh2particle(meshpos, config, meshdata, pos::AbstractVector{T}, ::VertexMode, ::TSC) where T<:Number
    id, r = TSCinfo(meshpos, config, pos)
    return _gather(meshdata, id, r, Val(config.dim), Val(3))
end

@inline function _gather(meshdata, id, r, ::Val{1}, ::Val{n}) where n
    s = zero(eltype(meshdata))
    @inbounds for i in 1:n
        s += meshdata[id[i][1]] * r[i][1]
    end
    return s
end

@inline function _gather(meshdata, id, r, ::Val{2}, ::Val{n}) where n
    s = zero(eltype(meshdata))
    @inbounds for j in 1:n, i in 1:n
        s += meshdata[id[i][1], id[j][2]] * r[i][1] * r[j][2]
    end
    return s
end

@inline function _gather(meshdata, id, r, ::Val{3}, ::Val{n}) where n
    s = zero(eltype(meshdata))
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        s += meshdata[id[i][1], id[j][2], id[k][3]] * r[i][1] * r[j][2] * r[k][3]
    end
    return s
end

mesh2particle(meshpos, config, meshdata, pos::AbstractPoint, args...) = mesh2particle(meshpos, config, meshdata, SVector(pos), args...)
mesh2particle(meshpos, config, meshdata, pos::AbstractPoint) = mesh2particle(meshpos, config, meshdata, SVector(pos), config.mode, config.assignment)

"""
$(TYPEDSIGNATURES)

Assign particle data by inverse mass assignment.

## Examples
```julia
assignparticle(data, m, :Acc, :acc)
```
"""
Base.@propagate_inbounds function assignparticle(particles::StructArray, mesh::MeshCartesianStatic, symbolParticle::Symbol, symbolMesh::Symbol)
    config = mesh.config
    meshfield = getproperty(mesh, symbolMesh)

    CUDA.@allowscalar meshpos = Array(mesh.pos)
    CUDA.@allowscalar meshdata = Array(meshfield.data)
    CUDA.@allowscalar particlepos = Array(particles.Pos)
    CUDA.@allowscalar particledata = Array(getproperty(particles, symbolParticle))

    for i in eachindex(particles)
        if is_inbound(particlepos[i], config)
            pos = SVector(particlepos[i])
            particledata[i] = mesh2particle(meshpos, config, meshdata, pos, config.mode, config.assignment)
        end
    end

    if getproperty(particles, symbolParticle) isa CuArray
        getproperty(particles, symbolParticle) .= cu(particledata)
    else
        getproperty(particles, symbolParticle) .= particledata
    end
end
