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

COLORS = Dict( 
    :red     => "\x1b[31m",
    :green   => "\x1b[32m",
    :yellow  => "\x1b[33m",
    :blue    => "\x1b[34m",
    :magenta => "\x1b[35m",
    :cyan    => "\x1b[36m",
    :white   => "\x1b[37m",
    :bold    => "\x1b[1m",
    :normal  => "\x1b[0m",
    )

# Alias to print with color
pcolor = print_with_color
#function pcolor(col::Symbol, msg::AbstractString...)
    #print(COLORS[col], msg..., COLORS[:normal])
#end

# Print bold with color
function pbcolor(col::Symbol, msg::AbstractString...)
    const BOLD    = "\x1b[1m"
    const DEFAULT = "\x1b[0m"
    if is_linux() print(BOLD) end
    print_with_color(col, msg...)
    if is_linux() print(DEFAULT) end
end

using Base

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

#Use non-registered (jet) package FemMesh
#Pkg.installed("FemMesh") == nothing && Pkg.clone("https://github.com/RaulDurand/FemMesh")

try
    eval(:(using FemMesh))
catch err
    println(err)
    Pkg.clone("https://github.com/RaulDurand/FemMesh")
end

@reexport using FemMesh
import FemMesh.save # to be extended
import FemMesh.update! # to be extended

using JSON

# Tools
include("tools/linalg.jl")
include("tools/expr.jl")
include("tools/table.jl")

# Fem module
include("node.jl")
include("elem.jl")
include("face.jl")
include("bcs.jl")
include("track.jl")
include("domain.jl")
include("mec/solver.jl")
include("mec/mechanical.jl")
#include("seep/seep.jl")

end#module
