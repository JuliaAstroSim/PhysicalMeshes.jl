module PhysicalMeshes

using PrecompileTools
using Unitful, UnitfulAstro
using Distributed
using LinearAlgebra
using StaticArrays
#using Decimals
using StructArrays
using DocStringExtensions

using PhysicalParticles
using AstroSimBase

using CUDA
macro hascuda(expr)
    return has_cuda() ? :($(esc(expr))) : :(nothing)
end
@hascuda begin
    using CUDA
    CUDA.allowscalar(false)
end

import Core: Number
import Base: +, -, *, /, show, real, length, iterate, intersect, zero, iszero
import Unitful: Units, FloatTypes
#import Decimals: Decimal, decimal
import PhysicalParticles: PVector2D, PVector, area, volume, mass_center

export
    # Base
    +, -, *, /,
    show, real, length, iterate,

    # Traits
    PositivelyOriented,
    NegativelyOriented,
    UnOriented,

    Fast, Exact,

    Interior, Exterior,
    OnEdge, OnFace,

    CellMode, VertexMode,
    # Periodic, Dirichlet, Vacuum,

    NGP, CIC, TSC,

    # Core
    orient,
    orientation,
    centroid, #center,
    midpoint,
    normal,

    circumcenter,
    #circumcenter_exact,

    incircle,
    #incircle_exact,
    insphere,
    #insphere_exact,

    interior,
    exterior,

    #decimal,
    #floatnumber,
    # Sphere
    Sphere,

    # Tools
    zoom, shrink,

    # Line
    AbstractLine, AbstractLine2D, AbstractLine3D,
    Line, Line2D,

    # Plane
    AbstractPlane, AbstractPlane3D,
    Plane,
    coplanar,

    # Polygon
    AbstractPolygon, AbstractPolygon2D, AbstractPolygon3D,
    Polygon2D, Polygon3D,
    polygon_rect, polygon_regular,
    isconvex,

    # Ray
    AbstractRay, AbstractRay2D, AbstractRay3D,
    Ray2D, Ray3D,
    intersect,
    reflect,
    
    # Cube
    AbstractCube, AbstractCube2D, AbstractCube3D,
    Cube, Cube2D,

    # Triangle
    AbstractTriangle,
    AbstractTriangle2D, AbstractTriangle3D,
    Triangle, Triangle2D,

    # Tetrahedron
    Tetrahedron,
    
    # Mesh
    particle2mesh, assignmesh,
    mesh2particle, assignparticle,
    is_inbound,
    outbound_list,
    total_mass,

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

abstract type AbstractPlane{T} <: AbstractGeometryType{T} end
abstract type AbstractPlane3D{T} <: AbstractPlane{T} end

abstract type AbstractPolygon{T} <: AbstractGeometryType{T} end
abstract type AbstractPolygon2D{T} <: AbstractPolygon{T} end
abstract type AbstractPolygon3D{T} <: AbstractPolygon{T} end

abstract type AbstractRay{T} <: AbstractGeometryType{T} end
abstract type AbstractRay2D{T} <: AbstractRay{T} end
abstract type AbstractRay3D{T} <: AbstractRay{T} end

abstract type AbstractCube{T} <: AbstractGeometryType{T} end
abstract type AbstractCube2D{T} <: AbstractCube{T} end
abstract type AbstractCube3D{T} <: AbstractCube{T} end

abstract type AbstractTriangle{T} <: AbstractGeometryType{T} end
abstract type AbstractTriangle2D{T} <: AbstractTriangle{T} end
abstract type AbstractTriangle3D{T} <: AbstractTriangle{T} end

abstract type AbstractTetrahedron{T} <: AbstractGeometryType{T} end

abstract type AbstractMesh{T} <: AbstractGeometryType{T} end
# abstract type AbstractMesh1D{T} <: AbstractMesh{T} end
# abstract type AbstractMesh2D{T} <: AbstractMesh{T} end
# abstract type AbstractMesh3D{T} <: AbstractMesh{T} end

include("Traits.jl")

#include("core/decimals.jl")
include("core/orient.jl")
include("core/circumcenter.jl")
include("core/predicates.jl")

include("Sphere.jl")

include("Line.jl")
include("Plane.jl")
include("Polygon.jl")
include("Ray.jl")
include("Cube.jl")
include("Triangle.jl")
include("Tetrahedron.jl")

include("mesh/Mesh.jl")
include("mesh/particle2mesh.jl")
include("mesh/tools.jl")

include("precompile.jl")

end