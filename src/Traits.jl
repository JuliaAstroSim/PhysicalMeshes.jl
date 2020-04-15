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