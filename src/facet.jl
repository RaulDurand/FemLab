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

abstract type Facet end

mutable struct Face<:Facet
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


mutable struct Edge<:Facet
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union{Element,Void}
    isedge::Bool
    function Edge(shape, nodes, ndim)
        this = new(shape, nodes, ndim)
        this.oelem  = nothing
        this.isedge = false
        return this
    end
end


function getindex{T<:Facet}(facets::Array{T,1}, cond::Expr)
    condm = fix_comparison_arrays(cond)
   
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Facet getindex: Invalid condition ", cond)
    end

    result = Array{T}(0)
    for facet in facets
        coords = nodes_coords(facet.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        if @static VERSION>v"0.6.0-rc1.0" ? Base.invokelatest(fun, x, y, z) : fun(x, y, z)
            push!(result, facet) 
        end
    end

    if length(result) == 0
        warn("No facets found that match: $cond\n")
    end

    return result
end

#getindex{T<:Facet}(facets::Array{T,1}, cond::AbstractString) = getindex(facets, parse(cond))


# Get all nodes from a collection of facets
function get_nodes{T<:Facet}(facets::Array{T,1})
    nodes = Set{Node}()
    for facet in facets
        for node in facet.nodes
            push!(nodes, node)
        end
    end
    return [node for node in nodes]
end

# Index operator for a collection of facets
function getindex{T<:Facet}(facets::Array{T,1}, s::Symbol)
    if s == :nodes
        return get_nodes(facets)
    end
    error("Facet getindex: Invalid symbol $s")
end
