using Test
using Unitful, UnitfulAstro
import Interpolations: LinearInterpolation

using PhysicalParticles
using PhysicalMeshes
using AstroSimBase
import PhysicalMeshes: Unitless2D, Unitless3D, Physical2D, Physical3D
import PhysicalParticles: Ball, assign_particles

using Distributed
using FFTW

include("testSphere.jl")
include("testLine.jl")
include("testTriangle.jl")
include("testCube.jl")
include("testMesh.jl")
include("testTetrahedron.jl")
include("testPolygon.jl")
include("testRay.jl")
include("testPlane.jl")
include("testCore.jl")
include("testField.jl")
include("testParallel.jl")
include("testMeshOps.jl")
# include("testPIC.jl")