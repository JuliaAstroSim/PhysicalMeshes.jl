function NGPinfo(meshpos, config, pos::AbstractVector{T}) where T<:Number
    id = round(Int, (pos .- config.Min) ./ config.Δ) .+ 1 .+ config.NG
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
    config = mesh.config
    idl = floor.(Int, (pos .- config.Min) ./ config.Δ .- 0.5) .+ 1 .+ config.NG
    idm = idl .+ 1
    idr = idm .+ 1
    
    pl = SVector(meshpos[idl...])
    
    rl = (pl .- pos .+ 1.5 * config.Δ) .^2 ./ 2 ./ config.Δ ./ config.Δ
    rr = (pl .- pos .+ 0.5 * config.Δ) .^2 ./ 2 ./ config.Δ ./ config.Δ
    rm = 1.0 .- rl .- rr    # assignment ratio of middle vertex

    id = [idl, idm, idr]
    r = [rl, rm, rr]
    return id, r
end

function is_inbound(pos::PVector, config::MeshConfig)
    if pos.x > config.Max[1] || pos.x < config.Min[1] ||
        pos.y > config.Max[2] || pos.y < config.Min[2] ||
        pos.z > config.Max[3] || pos.z < config.Min[3]
        return false
    else
        return true
    end
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

function particle2mesh!(meshpos, config, pos::AbstractVector{T}, meshdata, particledata, ::VertexMode, ::NGP) where T<:Number
    # Find the nearest vertex and assign with particledata
    id = NGPinfo(meshpos, config, pos)
    meshdata[id...] += particledata
end

function particle2mesh!(meshpos, config, pos::AbstractVector{T}, meshdata, particledata, ::VertexMode, ::CIC) where T<:Number
    id, r = CICinfo(meshpos, config, pos)
    for i in 1:2, j in 1:2, k in 1:2
        @inbounds meshdata[id[i][1], id[j][2], id[k][3]] += particledata * r[i][1] * r[j][2] * r[k][3]
    end
end

function particle2mesh!(meshpos, config, pos::AbstractVector{T}, meshdata, particledata, ::VertexMode, ::TSC) where T<:Number
    id, r = TSCinfo(meshpos, config, pos)
    for i in 1:3, j in 1:3, k in 1:3
        @inbounds meshdata[id[i][1], id[j][2], id[k][3]] += particledata * r[i][1] * r[j][2] * r[k][3]
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
    getproperty(mesh, symbolMesh) .*= 0.0

    CUDA.@allowscalar meshpos = Array(mesh.pos)
    CUDA.@allowscalar meshdata = Array(getproperty(mesh, symbolMesh))
    CUDA.@allowscalar particlepos = Array(particles.Pos)
    CUDA.@allowscalar particledata = Array(getproperty(particles, symbolParticle))

    for i in eachindex(particles)
        rho = particledata[i] / prod(config.Δ)
        if is_inbound(particlepos[i], config)
            pos = SVector(particlepos[i])
            particle2mesh!(meshpos, config, pos, meshdata, rho, config.mode, config.assignment)
        end
    end

    if getproperty(mesh, symbolMesh) isa CuArray
        getproperty(mesh, symbolMesh) .= cu(meshdata)
    else
        getproperty(mesh, symbolMesh) .= meshdata
    end
end

assignmesh(m::MeshCartesianStatic) = assignmesh(m.data, m, :Mass, :rho)

function mesh2particle(meshpos, config, meshdata, pos::AbstractVector{T}, ::VertexMode, ::NGP) where T<:Number
    id = NGPinfo(meshpos, config, pos)
    return meshdata[id...]
end

function mesh2particle(meshpos, config, meshdata, pos::AbstractVector{T}, ::VertexMode, ::CIC) where T<:Number
    id, r = CICinfo(meshpos, config, pos)
    return sum([meshdata[id[i][1], id[j][2], id[k][3]] * r[i][1] * r[j][2] * r[k][3] for i in 1:2, j in 1:2, k in 1:2])
end

function mesh2particle(meshpos, config, meshdata, pos::AbstractVector{T}, ::VertexMode, ::TSC) where T<:Number
    id, r = TSCinfo(meshpos, config, pos)
    return sum([meshdata[id[i][1], id[j][2], id[k][3]] * r[i][1] * r[j][2] * r[k][3] for i in 1:3, j in 1:3, k in 1:3])
end

mesh2particle(meshpos, config, meshdata, pos::AbstractPoint, args...) = mesh2particle(meshpos, config, meshdata, SVector(pos), args...)

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

    CUDA.@allowscalar meshpos = Array(mesh.pos)
    CUDA.@allowscalar meshdata = Array(getproperty(mesh, symbolMesh))
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
