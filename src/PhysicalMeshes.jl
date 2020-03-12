module PhysicalMeshes

using Unitful, UnitfulAstro
using Distributed

using PhysicalParticles

import Base: +, -, *, /, show, real, length, iterate
import Unitful: Units

export
    # Base
    +, -, *, /,
    show, real, length, iterate,

    # Line
    Line, Line2D,
    
    # Cube
    Cube, Cube2D,

    # Triangle
    Triangle, Triangle2D,

    # Tetrahedron

    
    # Mesh
    orientation,

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
include("Line.jl")
include("Cube.jl")
include("Triangle.jl")
include("Tetrahedron.jl")

include("PrettyPrinting.jl")

end