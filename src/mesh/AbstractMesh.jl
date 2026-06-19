abstract type AbstractMesh{T} <: AbstractGeometryType{T} end
abstract type AbstractMesh1D{T} <: AbstractMesh{T} end
abstract type AbstractMesh2D{T} <: AbstractMesh{T} end
abstract type AbstractMesh3D{T} <: AbstractMesh{T} end

# Abstract methods for all meshes
@inline function Base.getproperty(mesh::AbstractMesh, name::Symbol)
    if name in fieldnames(typeof(mesh))
        return Base.getfield(mesh, name)
    else
        error("Field $name not found in mesh")
    end
end

"""
$(TYPEDSIGNATURES)
Default `setproperty!` for abstract meshes. It uses the unexported
`Base.setfield!` directly to avoid recursing back into the
user-overridden `setproperty!` method (which would cause an infinite
loop on subtypes such as `MeshUnstructured` that want to dispatch
custom field storage).
"""
@inline function Base.setproperty!(mesh::AbstractMesh, name::Symbol, value)
    if name in fieldnames(typeof(mesh))
        return Core.setfield!(mesh, name, value)
    else
        error("Field $name not found in mesh")
    end
end

function Base.show(io::IO, mesh::AbstractMesh)
    print(io, typeof(mesh))
end

# Mesh operations
function cell_volume(mesh::AbstractMesh, cell_id::Int)
    error("cell_volume not implemented for $(typeof(mesh))")
end

function cell_area(mesh::AbstractMesh, cell_id::Int)
    error("cell_area not implemented for $(typeof(mesh))")
end

function cell_center(mesh::AbstractMesh, cell_id::Int)
    error("cell_center not implemented for $(typeof(mesh))")
end

function interpolate_field(mesh::AbstractMesh, field::AbstractField, pos::PVector)
    error("interpolate_field not implemented for $(typeof(mesh))")
end

function add_node!(mesh::AbstractMesh, node::PVector)
    error("add_node! not implemented for $(typeof(mesh))")
end

function add_cell!(mesh::AbstractMesh, cell::AbstractArray{Int})
    error("add_cell! not implemented for $(typeof(mesh))")
end

function add_field!(mesh::AbstractMesh, name::Symbol, field::AbstractField)
    error("add_field! not implemented for $(typeof(mesh))")
end

function get_field(mesh::AbstractMesh, name::Symbol)
    error("get_field not implemented for $(typeof(mesh))")
end

# Dimension-specific abstract methods
function node_neighbors(mesh::AbstractMesh2D, node_id::Int)
    error("node_neighbors not implemented for $(typeof(mesh))")
end

function node_neighbors(mesh::AbstractMesh3D, node_id::Int)
    error("node_neighbors not implemented for $(typeof(mesh))")
end

function cell_neighbors(mesh::AbstractMesh2D, cell_id::Int)
    error("cell_neighbors not implemented for $(typeof(mesh))")
end

function cell_neighbors(mesh::AbstractMesh3D, cell_id::Int)
    error("cell_neighbors not implemented for $(typeof(mesh))")
end

# Mesh utilities
function mesh_dimension(mesh::AbstractMesh1D)
    return 1
end

function mesh_dimension(mesh::AbstractMesh2D)
    return 2
end

function mesh_dimension(mesh::AbstractMesh3D)
    return 3
end

function mesh_dimension(mesh::AbstractMesh)
    return mesh.config.dim
end

function number_of_nodes(mesh::AbstractMesh)
    error("number_of_nodes not implemented for $(typeof(mesh))")
end

function number_of_cells(mesh::AbstractMesh)
    error("number_of_cells not implemented for $(typeof(mesh))")
end

function number_of_fields(mesh::AbstractMesh)
    error("number_of_fields not implemented for $(typeof(mesh))")
end

# Concrete-mesh implementations are added next to each struct definition to
# avoid forward-reference errors when the include order is re-arranged.
