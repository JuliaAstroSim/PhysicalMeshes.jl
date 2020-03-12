using Test
using Unitful, UnitfulAstro

using PhysicalParticles
using PhysicalMeshes
import PhysicalMeshes: Unitless2D, Unitless3D, Physical2D, Physical3D


include("testLine.jl")
include("testCube.jl")
include("testTriangle.jl")
include("testMesh.jl")
include("testTetrahedron.jl")