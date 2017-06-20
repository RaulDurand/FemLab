
include("mechanical.jl")
export elem_config_dofs, elem_init, elem_stiffness, elem_map, elem_RHS, elem_dF!, elem_and_node_vals
export set_state

# Models for solid elements (2D and 3D)
include("abs-solid.jl")
include("elasticsolid.jl")
include("dp.jl")
include("kotsovos.jl")
include("mazars.jl")
export AbsSolid, ElasticSolid, DruckerPrager, Kotsovos, Mazars

# Models for truss elements
include("abs-truss.jl")
include("truss.jl")
include("pptruss.jl")
export AbsTruss, Truss, PPTruss

# Models for embedded truss
include("abs-embtruss.jl")
include("embtruss.jl")
include("embpptruss.jl")
export AbsEmbTruss, EmbTruss, EmbPPTruss

# Models for joint elements
include("abs-joint.jl")
include("joint.jl")
include("mcjoint.jl")
export AbsJoint, Joint, MCJoint

# Models for 1D joint elements
include("abs-joint1d.jl")
include("joint1d.jl")
include("mcjoint1d.jl")
include("cebjoint1d.jl")
export AbsJoint1D, Joint1D, MCJoint1D, CEBJoint1D

include("solver.jl")
export solve!, MecSolverData

include("abs-beam.jl")
include("beam.jl")
