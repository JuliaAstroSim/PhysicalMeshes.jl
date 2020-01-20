module PhysicalMeshes

using Unitful, UnitfulAstro
using Distributed

using PhysicalParticles

import Base: show, real, length, iterate

export
    # Base
    show, real, length, iterate,

    # Line
    AbstractLine, AbstractLine2D, AbstractLine3D,
    
    # Triangle
    AbstractTriangle, AbstractTriangle2D, AbstractTriangle3D,

    # Tetrahedron

    
    # Mesh
    AbstractMesh, AbstractMesh2D, AbstractMesh3D,
    
    orientation,
    volume,
    area,
    len # Circumference


abstract type AbstractLine{T} end
abstract type AbstractLine2D{T} <: AbstractLine{T} end
abstract type AbstractLine3D{T} <: AbstractLine{T} end

@inline real(p::T) where T <: AbstractLine = p


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
include("Triangle.jl")
include("Tetrahedron.jl")

include("Voronoi/Voronoi.jl")

include("PrettyPrinting.jl")

end