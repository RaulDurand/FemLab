##############################################################################
#    FemLab - Finite Element Library                                         #
#    Copyright (C) 2014 Raul Durand <raul.durand at gmail.com>               #
#                                                                            #
#    This file is part of FemLab.                                            #
#                                                                            #
#    FemLab is free software: you can redistribute it and/or modify          #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    any later version.                                                      #
#                                                                            #
#    FemLab is distributed in the hope that it will be useful,               #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with FemLab.  If not, see <http://www.gnu.org/licenses/>.         #
##############################################################################

VERSION >= v"0.4.0-dev+6521" && __precompile__()

module FemLab

using Base

# Definitions module
include("definitions.jl")
include("tools/linalg.jl")
include("tools/table.jl")

export DTable , push!, getindex, save, loadtable

# Mesh module
include("mesh/mesh.jl")
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

#include("seep/seep.jl")

end#module
