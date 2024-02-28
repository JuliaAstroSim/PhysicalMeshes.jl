using Test
using Unitful, UnitfulAstro
import Interpolations: LinearInterpolation

using PhysicalParticles
using PhysicalMeshes
import PhysicalMeshes: Unitless2D, Unitless3D, Physical2D, Physical3D


include("testSphere.jl")
include("testLine.jl")
include("testTriangle.jl")
include("testCube.jl")
include("testMesh.jl")
include("testTetrahedron.jl")