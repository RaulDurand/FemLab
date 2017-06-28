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

__precompile__() 

"""
**FemLab.jl**

FemLab module implements functions and types to perform finite element analyses.

**Important data types**

Node, Element, Domain, Dof, Ip, NodeBC, FaceBC

**Important functions** 

set_mat, set_bc, clear_bc, solve!, save

"""
module FemLab
using  Base, JSON, Reexport

#Using non-registered (jet) package FemMesh
#Pkg.installed("FemMesh") == nothing && Pkg.clone("https://github.com/RaulDurand/FemMesh")

try
    eval(:(using FemMesh))
catch err
    #println(err)
    Pkg.clone("https://github.com/RaulDurand/FemMesh")
end

@reexport using FemMesh
import FemMesh.save # to be extended
import FemMesh.update! # to be extended

# Tools module
include("tools/linalg.jl")
include("tools/expr.jl")
include("tools/table.jl")
include("tools/tensors.jl")

# generic exports
export max, min, sort, reset, getindex, sort, copy!, show

# Fem module
include("globals.jl")

include("node.jl")
export Node, Dof, add_dof

include("ip.jl")
export Ip, ip_vals, maximum, minimum, sort

include("material.jl")
export Material, read_prms

include("element.jl")
export Element
export set_mat, get_nodes, get_ips, set_state, reset, getcoords, getvals, get_map

include("facet.jl")
export Facet, Face, Edge

include("bcs.jl")
export NodeBC, FaceBC, EdgeBC, apply_bc

include("monitors.jl")
export NodeMonitor, NodesMonitor, IpMonitor, IpsMonitor, FacesMonitor, EdgesMonitor
export NodeTracker, NodesTracker, IpTracker, IpsTracker, FacesTracker, EdgesTracker
export update_monitor!

include("domain.jl")
export set_bc, set_monitors, update_monitors!, set_trackers

include("io.jl")
export show


# Mechanical module
include("mech/include.jl")


# Porous media flow module
#include("seep/include.jl")

end#module
