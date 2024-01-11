abstract type AbstractOrientation end
struct PositivelyOriented <: AbstractOrientation end
struct NegativelyOriented <: AbstractOrientation end
struct UnOriented <: AbstractOrientation end

abstract type PrecisionMode end
struct Fast <: PrecisionMode end
struct Exact <: PrecisionMode end

abstract type AbstractPredicate end
struct Interior <: AbstractPredicate end
struct Exterior <: AbstractPredicate end
struct OnEdge <: AbstractPredicate end
struct OnFace <: AbstractPredicate end

abstract type MeshMode end
struct CellMode <: MeshMode end
struct VertexMode <: MeshMode end

@inline length(p::T) where T <: MeshMode = 1
@inline iterate(p::T) where T <: MeshMode = (p,nothing)
@inline iterate(p::T,st) where T <: MeshMode = nothing

abstract type MeshAssignment end
struct NGP <: MeshAssignment end # nearest grid point
struct CIC <: MeshAssignment end # cloud in cell
struct TSC <: MeshAssignment end # triangular shaped cloud

@inline length(p::T) where T <: MeshAssignment = 1
@inline iterate(p::T) where T <: MeshAssignment = (p,nothing)
@inline iterate(p::T,st) where T <: MeshAssignment = nothing

