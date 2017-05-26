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

type Face
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union{Element,Void}
    isedge::Bool
    function Face(shape, nodes, ndim)
        this = new(shape, nodes, ndim)
        this.oelem  = nothing
        this.isedge = false
        return this
    end
end

const Edge=Face

function getindex(faces::Array{Face,1}, cond::Expr)
    condm = fix_comparison_arrays(cond)
   
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Faces getindex: Invalid condition ", cond)
    end

    result = Array{Face}(0)
    for face in faces
        coords = getcoords(face.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        if @static VERSION>v"0.6.0-rc1.0" ? Base.invokelatest(fun, x, y, z) : fun(x, y, z)
            push!(result, face) 
        end
    end

    if length(result) == 0
        printcolor(:red, "Warning: No faces found that match: $cond\n")
    end

    return result
end

getindex(faces::Array{Face,1}, cond::AbstractString) = getindex(faces, parse(cond))


# Get all nodes from a collection of faces
function get_nodes(faces::Array{Face,1})
    nodes = Set{Node}()
    for face in faces
        for node in face.nodes
            push!(nodes, node)
        end
    end
    return [node for node in nodes]
end

# Index operator for a collection of faces
function getindex(faces::Array{Face,1}, s::Symbol)
    if s == :nodes
        return get_nodes(faces)
    end
    error("Face getindex: Invalid symbol $s")
end
