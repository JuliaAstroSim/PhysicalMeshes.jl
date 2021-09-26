module PhysicalMeshes

using Reexport
using Unitful, UnitfulAstro
using Distributed
using LinearAlgebra
using StaticArrays
#using Decimals
using StructArrays

@reexport using PhysicalParticles

import Core: Number
import Base: +, -, *, /, show, real, length, iterate
import Unitful: Units, FloatTypes
#import Decimals: Decimal, decimal
import PhysicalParticles: PVector2D, PVector, area, volume

export
    # Base
    +, -, *, /,
    show, real, length, iterate,

    # Traits
    PositivelyOriented,
    NegativelyOriented,
    UnOriented,

    Fast, Exact,

    Inner, Outer,
    OnEdge, OnFace,

    CellMode, VertexMode,

    NGP, CIC, TSC,

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
    #circumcenter_exact,

    incircle,
    #incircle_exact,
    insphere,
    #insphere_exact,

    #decimal,
    #floatnumber,

    # Line
    AbstractLine,
    AbstractLine2D, AbstractLine3D,
    Line, Line2D,
    
    # Cube
    AbstractCube,
    AbstractCube2D, AbstractCube3D,
    Cube, Cube2D,

    # Triangle
    AbstractTriangle,
    AbstractTriangle2D, AbstractTriangle3D,
    Triangle, Triangle2D,

    # Tetrahedron
    Tetrahedron,
    
    # Mesh
    centroid, #center,
    midpoint,

    AbstractMesh,
    AbstractMesh3D, AbstractMesh2D, AbstractMesh1D,
    MeshConfig,
    MeshCartesianStatic,

    #volume,
    #area,
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
abstract type AbstractMesh1D{T} <: AbstractMesh{T} end
abstract type AbstractMesh2D{T} <: AbstractMesh{T} end
abstract type AbstractMesh3D{T} <: AbstractMesh{T} end

include("Traits.jl")

#include("core/decimals.jl")
include("core/orient.jl")
include("core/circumcenter.jl")
include("core/predicates.jl")

include("Line.jl")
include("Cube.jl")
include("Triangle.jl")
include("Tetrahedron.jl")

include("mesh/Mesh.jl")
include("mesh/particle2mesh.jl")

include("PrettyPrinting.jl")

end