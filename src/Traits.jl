abstract type AbstractMeshType end
struct Physical2D <: AbstractMeshType end
struct Physical3D <: AbstractMeshType end
struct Unitless2D <: AbstractMeshType end
struct Unitless3D <: AbstractMeshType end

meshtype(p::AbstractParticle) = meshtype(p.Pos, p.Pos.x)
meshtype(p::AbstractPoint) = meshtype(p, p.x)

meshtype(::AbstractPoint2D, ::Number) = Unitless2D()
meshtype(::AbstractPoint3D, ::Number) = Unitless3D()
meshtype(::AbstractPoint2D, ::Quantity) = Physical2D()
meshtype(::AbstractPoint3D, ::Quantity) = Physical3D()

meshtype(a::Array) = meshtype(a[1])

abstract type AbstractOrientation end
struct PositivelyOriented <: AbstractOrientation end
struct NegativelyOriented <: AbstractOrientation end