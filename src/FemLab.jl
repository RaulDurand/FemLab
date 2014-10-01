module FemLab

using Base

# Definitions module
include("definitions.jl")
using  .Definitions
export DTable , push!, getindex, save, loadtable

# Mesh module
include("mesh/meshgen.jl")
#using .MeshGen
export is_solid, is_line, is_joint
export Block2D, Block3D, BlockTruss, BlockInset, Mesh
export generate_mesh#, save

# Fem module
include("node.jl")
include("elem.jl")
include("domain.jl")
include("mec/solver.jl")
include("mec/mechanical.jl")


end#module
