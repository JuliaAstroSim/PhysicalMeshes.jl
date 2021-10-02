function particle2mesh(mesh::AbstractMesh, pos::AbstractArray, rho::Number, ::VertexMode, ::NGP)
    
end

function particle2mesh(mesh::AbstractMesh, pos::AbstractArray, rho::Number, ::VertexMode, ::CIC)
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

    for i in 1:2, j in 1:2, k in 1:2     # 
        mesh.rho[id[i][1], id[j][2], id[k][3]] = rho * r[i][1] * r[j][2] * r[k][3]
    end
end

function particle2mesh(mesh::AbstractMesh, ::TSC)
    
end

function assignmesh(particles::StructArray, mesh::MeshCartesianStatic, ::NGP)
    
end

function assignmesh(particles::StructArray, mesh::MeshCartesianStatic, ::CIC)
    config = mesh.config
    @batch for i in eachindex(particles)
        rho = particles.Mass[i] / prod(config.Δ)
        pos = SVector(particles.Pos[i])
        particle2mesh(mesh, pos, rho, config.mode, config.assignment)
    end
end

function mesh2particle()
    
end