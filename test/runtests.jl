using Test
using Unitful, UnitfulAstro

using PhysicalParticles
using PhysicalMeshes
import PhysicalMeshes: Unitless2D, Unitless3D, Physical2D, Physical3D


include("testLine.jl")
include("testTriangle.jl")
include("testCube.jl")
include("testMesh.jl")
include("testTetrahedron.jl")

#=
include("PhysicalMeshes.jl\\src\\PhysicalMeshes.jl"); using Main.PhysicalMeshes, Unitful, PhysicalParticles

=#