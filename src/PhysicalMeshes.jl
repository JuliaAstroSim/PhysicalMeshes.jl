module PhysicalMeshes

using Unitful, UnitfulAstro
using Distributed

using PhysicalParticles

import Base: show, real, length, iterate

export
    # Base
    show, real, length, iterate,

    # Line
    Line, Line2D,
    
    # Cube
    Cube, Cube2D,

    # Triangle
    Triangle, Triangle2D,

    # Tetrahedron

    
    # Mesh
    AbstractMesh, AbstractMesh2D, AbstractMesh3D,
    
    orientation,
    volume,
    area,
    centroid,
    len # Circumference


abstract type AbstractLine{T} end
abstract type AbstractLine2D{T} <: AbstractLine{T} end
abstract type AbstractLine3D{T} <: AbstractLine{T} end

@inline real(p::T) where T <: AbstractLine = p

abstract type AbstractCube{T} end
abstract type AbstractCube2D{T} <: AbstractCube{T} end
abstract type AbstractCube3D{T} <: AbstractCube{T} end

@inline real(p::T) where T <: AbstractCube = p

abstract type AbstractTriangle{T} end
abstract type AbstractTriangle2D{T} <: AbstractTriangle{T} end
abstract type AbstractTriangle3D{T} <: AbstractTriangle{T} end

@inline real(p::T) where T <: AbstractTriangle = p


abstract type AbstractTetrahedron{T} end

@inline real(p::T) where T <: AbstractTetrahedron = p


abstract type AbstractMesh{T} end
abstract type AbstractMesh2D{T} <: AbstractMesh{T} end
abstract type AbstractMesh3D{T} <: AbstractMesh{T} end

@inline real(p::T) where T <: AbstractMesh = p


include("Traits.jl")
include("Line.jl")
include("Cube.jl")
include("Triangle.jl")
include("Tetrahedron.jl")

include("FVM/Voronoi.jl")

include("PrettyPrinting.jl")

end