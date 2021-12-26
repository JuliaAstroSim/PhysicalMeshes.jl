function NGPinfo(mesh::AbstractMesh, pos::AbstractVector{T}) where T<:Number
    config = mesh.config
    id = round(Int, (pos .- config.Min) ./ config.Δ) .+ 1 .+ config.NG
end

function CICinfo(mesh::AbstractMesh, pos::AbstractVector{T}) where T<:Number
    config = mesh.config
    # Here, left means smaller coordinate value in that dimension
    # find the index of nearest left and right vertex
    idl = floor.(Int, (pos .- config.Min) ./ config.Δ) .+ 1 .+ config.NG
    idr = idl .+ 1

    pl = SVector(mesh.pos[idl...])       # position of left vertex

    rl = (pl .- pos) ./ config.Δ .+ 1.0  # assignment ratio of left vertex
    rr = 1.0 .- rl                       # assignment ratio of right vertex

    id = [idl, idr]                      # pack info into array
    r = [rl, rr]

    return id, r
end

function TSCinfo(mesh::AbstractMesh, pos::AbstractVector{T}) where T<:Number
    config = mesh.config
    idl = floor.(Int, (pos .- config.Min) ./ config.Δ .- 0.5) .+ 1 .+ config.NG
    idm = idl .+ 1
    idr = idm .+ 1
    
    pl = SVector(mesh.pos[idl...])
    
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

outbound_list(m::MeshCartesianStatic) = outbound_list(m.data.Pos, m)

function particle2mesh(mesh::AbstractMesh, pos::AbstractVector{T}, rho::Number, symbolMesh::Symbol, ::VertexMode, ::NGP) where T<:Number
    # Find the nearest vertex and assign with rho
    id = NGPinfo(mesh, pos)
    getproperty(mesh, symbolMesh)[id...] += rho
end

function particle2mesh(mesh::AbstractMesh, pos::AbstractVector{T}, rho::Number, symbolMesh::Symbol, ::VertexMode, ::CIC) where T<:Number
    id, r = CICinfo(mesh, pos)
    for i in 1:2, j in 1:2, k in 1:2
        @inbounds getproperty(mesh, symbolMesh)[id[i][1], id[j][2], id[k][3]] += rho * r[i][1] * r[j][2] * r[k][3]
    end
end

function particle2mesh(mesh::AbstractMesh, pos::AbstractVector{T}, rho::Number, symbolMesh::Symbol, ::VertexMode, ::TSC) where T<:Number
    id, r = TSCinfo(mesh, pos)
    for i in 1:3, j in 1:3, k in 1:3
        @inbounds getproperty(mesh, symbolMesh)[id[i][1], id[j][2], id[k][3]] += rho * r[i][1] * r[j][2] * r[k][3]
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

    #TODO this is unsafe assign
    #Threads.@threads for k in 1:Threads.nthreads()
    #    Head, Tail = split_block(length(particles), k, Threads.nthreads())
    #    for i in Head:Tail
    #        rho = getproperty(particles, symbolParticle)[i] / prod(config.Δ)
    #        if is_inbound(particles.Pos[i], config)
    #            pos = SVector(particles.Pos[i])
    #            particle2mesh(mesh, pos, rho, symbolMesh, config.mode, config.assignment)
    #        end
    #    end
    #end

    for i in eachindex(particles)
        rho = getproperty(particles, symbolParticle)[i] / prod(config.Δ)
        if is_inbound(particles.Pos[i], config)
            pos = SVector(particles.Pos[i])
            particle2mesh(mesh, pos, rho, symbolMesh, config.mode, config.assignment)
        end
    end

    #TODO GPU
end

assignmesh(m::MeshCartesianStatic) = assignmesh(m.data, m, :Mass, :rho)

function mesh2particle(mesh::AbstractMesh, pos::AbstractVector{T}, symbolMesh::Symbol, ::VertexMode, ::NGP) where T<:Number
    id = NGPinfo(mesh, pos)
    return getproperty(mesh, symbolMesh)[id...]
end

function mesh2particle(mesh::AbstractMesh, pos::AbstractVector{T}, symbolMesh::Symbol, ::VertexMode, ::CIC) where T<:Number
    id, r = CICinfo(mesh, pos)
    return sum([getproperty(mesh, symbolMesh)[id[i][1], id[j][2], id[k][3]] * r[i][1] * r[j][2] * r[k][3] for i in 1:2, j in 1:2, k in 1:2])
end

function mesh2particle(mesh::AbstractMesh, pos::AbstractVector{T}, symbolMesh::Symbol, ::VertexMode, ::TSC) where T<:Number
    id, r = TSCinfo(mesh, pos)
    return sum([getproperty(mesh, symbolMesh)[id[i][1], id[j][2], id[k][3]] * r[i][1] * r[j][2] * r[k][3] for i in 1:3, j in 1:3, k in 1:3])
end

mesh2particle(mesh::AbstractMesh, pos::AbstractPoint, args...) = mesh2particle(mesh, SVector(pos), args...)

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
    Threads.@threads for k in 1:Threads.nthreads()
        Head, Tail = split_block(length(particles), k, Threads.nthreads())
        for i in Head:Tail
            if is_inbound(particles.Pos[i], config)
                pos = SVector(particles.Pos[i])
                getproperty(particles, symbolParticle)[i] = mesh2particle(mesh, pos, symbolMesh, config.mode, config.assignment)
            end
        end
    end
    #TODO GPU
end