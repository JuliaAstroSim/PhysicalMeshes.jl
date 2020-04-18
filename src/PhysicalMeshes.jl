module PhysicalMeshes

using Unitful, UnitfulAstro
using Distributed
using LinearAlgebra
using StaticArrays
using Decimals

using PhysicalParticles

import Base: +, -, *, /, show, real, length, iterate
import Unitful: Units, FloatTypes
import Decimals: decimal, number

export
    # Base
    +, -, *, /,
    show, real, length, iterate,

    # Traits
    PositivelyOriented,
    NegativelyOriented,
    UnOriented,

    Fast,
    Exact,

    Inner,
    Outter,
    OnEdge,
    OnFace,

    # Core
    orient,
    orientation,

    circumcenter,
    circumcenter_exact,

    incircle,
    incircle_exact,
    insphere,
    insphere_exact,

    decimal,
    number,

    # Line
    Line, Line2D,
    
    # Cube
    Cube, Cube2D,

    # Triangle
    Triangle, Triangle2D,

    # Tetrahedron
    Tetrahedron,
    
    # Mesh
    centroid, center, midpoint,

    volume,
    area,
    len # Circumference

abstract type AbstractGeometryType{T} end

@inline length(p::T) where T <: AbstractGeometryType = 1
@inline iterate(p::T) where T <: AbstractGeometryType = (p,nothing)
@inline iterate(p::T,st) where T <: AbstractGeometryType = nothing
@inline real(p::T) where T <: AbstractGeometryType = p

abstract type AbstractLine{T} <: AbstractGeometryType{T} end
abstract type AbstractLine2D{T} <: AbstractLine{T} end
abstract type AbstractLine3D{T} <: AbstractLine{T} end

abstract type AbstractCube{T} <: AbstractGeometryType{T} end
abstract type AbstractCube2D{T} <: AbstractCube{T} end
abstract type AbstractCube3D{T} <: AbstractCube{T} end

abstract type AbstractTriangle{T} <: AbstractGeometryType{T} end
abstract type AbstractTriangle2D{T} <: AbstractTriangle{T} end
abstract type AbstractTriangle3D{T} <: AbstractTriangle{T} end

abstract type AbstractTetrahedron{T} <: AbstractGeometryType{T} end

abstract type AbstractMesh{T} <: AbstractGeometryType{T} end
abstract type AbstractMesh2D{T} <: AbstractMesh{T} end
abstract type AbstractMesh3D{T} <: AbstractMesh{T} end

include("Traits.jl")

include("core/decimals.jl")
include("core/orient.jl")
include("core/circumcenter.jl")
include("core/predicates.jl")

include("Line.jl")
include("Cube.jl")
include("Triangle.jl")
include("Tetrahedron.jl")

include("PrettyPrinting.jl")

end