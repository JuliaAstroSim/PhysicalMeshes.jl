abstract type AbstractOrientation end
struct PositivelyOriented <: AbstractOrientation end
struct NegativelyOriented <: AbstractOrientation end
struct UnOriented <: AbstractOrientation end

abstract type PrecisionMode end
struct Fast <: PrecisionMode end
struct Exact <: PrecisionMode end

abstract type AbstractPredicate end
struct Inner <: AbstractPredicate end
struct Outter <: AbstractPredicate end
struct OnEdge <: AbstractPredicate end
struct OnFace <: AbstractPredicate end

abstract type MeshMode end
struct CellMode <: MeshMode end
struct VertexMode <: MeshMode end

abstract type MeshAssignment end
struct NGP <: MeshAssignment end # nearest grid point
struct CIC <: MeshAssignment end # cloud in cell
struct TSC <: MeshAssignment end # triangular shaped cloud