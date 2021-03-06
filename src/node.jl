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

import Base.reset
import Base.getindex
import Base.maximum
import Base.minimum
import Base.sort

# Dof
# ===

"""
`Dof()` 

Creates an object that represents a Degree of Freedom in a finite element analysis.
`Node` objects include a field called `dofs` which is an array of `Dof` objects.
"""
mutable struct Dof
    sU    ::Symbol  # essential bc name
    sF    ::Symbol  # natural bc name
    U     ::Float64 # essential value
    F     ::Float64 # natural value
    bryU  ::Float64
    bryF  ::Float64
    eq_id ::Int64   # number of equation in global system
    prescU::Bool    # flag for prescribed dof
    function Dof(sU::Symbol, sF::Symbol) 
        new(sU, sF, 0.0, 0.0, 0.0, 0.0, 0, false)
    end
end


# Node
# ====

"""
`Node(X)` 

Creates an object that represents a Node in a finite element analysis. The `X` parameter is a 
vector that represents the node coordinates.

**Important fields are**
`X`   : A vector of coordinates
`tag` : A string tag
`dofs`: An array of `Dof` objects
"""
mutable struct Node
    X       ::Array{Float64,1}
    tag     ::AbstractString
    id      ::Int
    dofs    ::Array{Dof,1}
    dofdict ::Dict{Symbol,Dof}
    n_shares::Int
 
    function Node(X::Array{Float64,1}; tag="", id=-1)
        this = new(X, tag, id)
        this.dofs = []
        this.dofdict = Dict()
        return this
    end
    function Node(point::Point; tag="", id=-1)
        Node([point.x, point.y, point.z], tag=tag, id=id)
    end
end

# Reset node information
function reset(nodes::Array{Node,1})
    for node in nodes
        for dof in node.dofs
            dof.U = 0.0
            dof.F = 0.0
            #dof.prescU = false
        end
    end
end


# Add a new degree of freedom to a node
function add_dof(node::Node, sU::Symbol, sF::Symbol)
    if !haskey(node.dofdict, sU)
        dof = Dof(sU, sF)
        push!(node.dofs, dof)
        node.dofdict[sU] = dof
        node.dofdict[sF] = dof
    end
end


# Index operator for node to get a dof
function getindex(node::Node, s::Symbol)
    return node.dofdict[s]
end


# Get node values in a dictionary
function node_vals(node::Node)
    coords = Dict( :x => node.X[1], :y => node.X[2], :z => node.X[3] )
    uvals  = Dict( dof.sU => dof.U for dof in node.dofs )
    fvals  = Dict( dof.sF => dof.F for dof in node.dofs )
    return merge(coords, uvals, fvals)
end


# Node collection
# ===============

# Index operator for an collection of nodes
function getindex(nodes::Array{Node,1}, cond::Expr) 
    condm = fix_comparison_scalar(cond)
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Node getindex: Invalid condition ", cond)
    end

    result = Array{Node}(0)
    for node in nodes
        #if @static VERSION>v"0.6.0-rc1.0" ? Base.invokelatest(fun, node.X[1], node.X[2], node.X[3]) : fun(node.X[1], node.X[2], node.X[3])
        if Base.invokelatest(fun, node.X[1], node.X[2], node.X[3])
            push!(result, node)
        end
    end

    if length(result) == 0
        warn("No nodes found that match: $cond\n")
    end

    return result
end

getindex(nodes::Array{Node,1}, cond::AbstractString) = getindex(nodes, parse(cond))


# Get node coordinates for an collection of nodes as a matrix
function nodes_coords(nodes::Array{Node,1}, ndim=3)
    nnodes = length(nodes)
    [ nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end


# Get the maximum value of a given coordinate for the whole collection of nodes
function maximum(nodes::Array{Node,1}, dir::Symbol=:x)
    idx = findfisrt((:x, :y, :z), dir)
    maximum([node.X[idx] for node in nodes])
end


function minimum(nodes::Array{Node,1}, dir::Symbol=:x)
    idx = findfisrt((:x, :y, :z), dir)
    minimum([node.X[idx] for node in nodes])
end


# Sort a collection of nodes in a given direction
function sort(nodes::Array{Node,1}, dir::Symbol=:x; rev::Bool=false)
    idx  = findfirst((:x, :y, :z), dir)
    idxs = sortperm([node.X[idx] for node in nodes], rev=rev)
    return nodes[idxs]
end

