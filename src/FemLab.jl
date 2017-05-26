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

#VERSION >= v"0.4.0-dev+6521" && __precompile__()
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
using  Base, JSON

# Code snippet from simonster
# https://github.com/simonster/Reexport.jl
macro reexport(ex)
    isa(ex, Expr) && (ex.head == :module ||
                      ex.head == :using ||
                      (ex.head == :toplevel &&
                       all(e->isa(e, Expr) && e.head == :using, ex.args))) ||
        error("@reexport: syntax error")

    if ex.head == :module
        modules = Any[ex.args[2]]
        ex = Expr(:toplevel, ex, Expr(:using, :., ex.args[2]))
    elseif ex.head == :using
        modules = Any[ex.args[end]]
    else
        modules = Any[e.args[end] for e in ex.args]
    end

    esc(Expr(:toplevel, ex,
             [:(eval(Expr(:export, setdiff(names($(mod)), [mod])...))) for mod in modules]...))
end

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

include("element.jl")
export Element, Ip
export set_mat, get_nodes, get_ips, set_state, reset, getcoords, getvals, get_map

include("face.jl")
export Face

include("bcs.jl")
export NodeBC, FaceBC, EdgeBC, apply_bc

include("monitor.jl")
export NodeMonitor, NodesMonitor, IpMonitor, IpsMonitor, FacesMonitor, EdgesMonitor
export NodeTracker, NodesTracker, IpTracker, IpsTracker, FacesTracker, EdgesTracker

include("domain.jl")
export set_bc, set_monitors, update_monitors, set_trackers


# Mechanical module
include("mech/include.jl")


# Porous media flow module
#include("seep/include.jl")

end#module
