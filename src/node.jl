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
import Base.max
import Base.min
import Base.sort

export Node
export set_bc, max, min, sort, reset


### Dof type

type Dof
    sU::Symbol
    sF::Symbol
    U ::Float64
    F ::Float64
    bryU::Float64
    bryF::Float64
    eq_id::Int64
    prescU::Bool
    function Dof(sU::Symbol, sF::Symbol)
        new(sU, sF, 0.0, 0.0, 0.0, 0.0, 0, false)
    end
end


### Node type


type Node
    X::Array{Float64,1}
    tag::String
    id ::Int
    dofs::Array{Dof,1}
    dofdict::Dict{Symbol,Dof}
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

function reset(nodes::Array{Node,1})
    for node in nodes
        for dof in node.dofs
            dof.U = 0.0
            dof.F = 0.0
        end
    end
end


function add_dof(node::Node, sU::Symbol, sF::Symbol)
    if !haskey(node.dofdict, sU)
        dof = Dof(sU, sF)
        push!(node.dofs, dof)
        node.dofdict[sU] = dof
        node.dofdict[sF] = dof
    end
end

function getindex(node::Node, s::Symbol)
    return node.dofdict[s]
end

function getvals(node::Node)
    uvals = [ dof.sU => dof.U for dof in node.dofs]
    fvals = [ dof.sF => dof.F for dof in node.dofs]
    merge(uvals, fvals)
end

function set_bc(node::Node; args...)
    for (key,val) in args
        if !haskey(node.dofdict, key); error("key ($key) not found in node ($(node.id)).") end
        dof = node.dofdict[key]
        if key==dof.sU
            dof.prescU = true
            dof.bryU  = val
        else
            dof.bryF += val
        end
    end
end

function clear_bc(node::Node)
    for dof in node.dofs
        dof.bryU   = 0.0
        dof.bryF   = 0.0
        dof.prescU = false
        dof.eq_id  = -1
    end
end


### Node collection


function getindex(nodes::Array{Node,1}, cond::Expr) 
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = cond
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Node getindex: Invalid condition ", cond)
    end

    result = Array(Node, 0)
    for node in nodes
        if fun(node.X[1], node.X[2], node.X[3])
            push!(result, node)
        end
    end
    return result
end

getindex(nodes::Array{Node,1}, cond::String) = getindex(nodes, parse(cond))


function getcoords(nodes::Array{Node,1}, ndim=3)
    nnodes = length(nodes)
    [ nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end

function set_bc(nodes::Array{Node,1}; args...) 
    if length(nodes)==0
        println("Warning, applying boundary conditions to empty array of nodes")
    end

    for node in nodes
        set_bc(node; args...)
    end
end

function max(nodes::Array{Node,1}, dir::Symbol) 
    idx = findfisrt((:x, :y, :z), dir)
    max([node.x[idx] for node in nodes]...)
end

function min(nodes::Array{Node,1}, dir::Symbol) 
    idx = findfisrt((:x, :y, :z), dir)
    min([node.x[idx] for node in nodes]...)
end

function sort(nodes::Array{Node,1}, dir::Symbol=:x, rev::Bool=false) 
    idx  = findfirst((:x, :y, :z), dir)
    idxs = sortperm([node.X[idx] for node in nodes], rev=rev)
    return nodes[idxs]
end

