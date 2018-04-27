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
import Base.copy!


@enum(ElemClass,
TRUSS_ELEM    = 1,
FRAME_ELEM    = 2,
SOLID_ELEM    = 3,
JOINT_ELEM    = 4,
JOINT1D_ELEM  = 5,
EMBEDDED_ELEM = 6
)

# Export
for s in instances(ElemClass)
    @eval export $(Symbol(s))
end 

"""
`Element(shape, nodes, ndim, [tag="",])`

Creates an 'Element' object for finite element analyses based on a
`shape`, an array of `nodes` and the space dimension `ndim`.
"""
mutable struct Element
    shape ::ShapeType
    nodes ::Array{Node,1}
    ndim  ::Int
    tag   ::TagType
    id    ::Int
    active::Bool
    ips   ::Array{Ip,1}
    mat   ::Material
    linked_elems::Array{Element,1}
    cache::Dict{Symbol,Any}

    function Element(shape, nodes, ndim, tag="")
        this        = new(shape, nodes, ndim, tag)
        this.active = true
        this.ips    = []
        this.linked_elems = []
        # Set the element class. For embedded elems, the class will be set by the Domain object.
        #this.class = shape.class
        this.cache = Dict{Symbol,Any}()
        this
    end
end

function reset(elems::Array{Element,1})
    # Resets data at integration points
    ndim = elems[1].ndim
    for elem in elems
        for ip in elem.ips
            ip.data = typeof(ip.data)(ndim)
        end
    end
end

# Get the element coordinates matrix
function elem_coords(elem::Element)
    nnodes = length(elem.nodes)
    ndim   = elem.ndim
    return [ elem.nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end



"""
`set_mat(elem, mat, [nips=0])`

Especifies the material model `mat` to be used to represent the behavior of an `Element` object `elem`.
"""
function set_mat(elem::Element, mat::Material; nips::Int64=0)
    # check if material is suitable for the element TODO: Improve this!
    embedded = elem.shape.class==LINE_SHAPE && length(elem.linked_elems)>0
    if embedded && client_shape_class(mat) != EMBEDDED
        error("set_mat: Material $(typeof(mat)) is not suitable for embedded elements")
    end

    # check if material is suitable for the element TODO: Improve this!
    if client_shape_class(mat) != elem.shape.class && !embedded
        error("set_mat: Material $(typeof(mat)) is not suitable for element class $(elem.shape.class)")
    end

    ipc =  get_ip_coords(elem.shape, nips)
    nips = size(ipc,1)

    resize!(elem.ips, nips)
    elem.mat = mat
    for i=1:nips
        R = ipc[i,1:3]
        w = ipc[i,4]
        elem.ips[i] = Ip(R, w)
        elem.ips[i].id = i
        elem.ips[i].data = new_ip_state(mat, elem.ndim)
        elem.ips[i].owner = elem
    end

    # finding ips global coordinates
    C     = elem_coords(elem)
    shape = elem.shape

    # fix for link elements
    if shape.class==JOINT1D_SHAPE
        bar   = elem.linked_elems[2]
        C     = elem_coords(bar)
        shape = bar.shape
    end

    # fix for joint elements
    if shape.class==JOINT_SHAPE
        C     = C[1:div(end,2),:]
        shape = shape.facet_shape
    end 

    # interpolation
    for ip in elem.ips
        N = shape.func(ip.R)
        ip.X = C'*N
        if length(ip.X)==2 # complete z=0.0 for 2D analyses
            ip.X = [ ip.X; 0.0 ]
        end
    end

    # configure degrees of freedom
    elem_config_dofs(mat, elem)

    # initialize element
    elem_init(mat, elem)
end



"""
`set_mat(elems, mat, [nips=0])`

Especifies the material model `mat` to be used to represent the behavior of a set of `Element` objects `elems`.
"""
function set_mat(elems::Array{Element,1}, mat::Material; nips::Int64=0)
    if length(elems)==0
        warn("Defining material model $(typeof(mat)) for an empty array of elements.\n")
    end
    for elem in elems
        set_mat(elem, mat, nips=nips)
    end
end


# Define the state at all integration points in a collection of elements
function set_state(elems::Array{Element,1}; args...)
    for elem in elems
        for ip in elem.ips
            set_state(ip.data; args...)
        end
    end
end


# Get all nodes from a collection of elements
function get_nodes(elems::Array{Element,1})
    nodes = Set{Node}()
    for elem in elems
        for node in elem.nodes
            push!(nodes, node)
        end
    end
    return [node for node in nodes]
end


# Get all ips from a collection of elements
function get_ips(elems::Array{Element,1})
    ips = Ip[]
    for elem in elems
        for ip in elem.ips
            push!(ips, ip)
        end
    end
    return ips
end

# Index operator for an element
function getindex(elem::Element, s::Symbol)
    if s == :nodes
        return elem.nodes
    end
    if s == :ips
        return elem.ips
    end
    error("Element getindex: Invalid symbol $s")
end

# Index operator for a collection of elements
function getindex(elems::Array{Element,1}, s::Symbol)
    if s == :solids
        return filter(elem -> elem.shape.class==SOLID_SHAPE, elems)
    end
    if s == :lines
        return filter(elem -> elem.shape.class==LINE_SHAPE, elems)
    end
    if s == :embedded
        return filter(elem -> elem.shape.class==LINE_SHAPE && length(elem.linked_elems)>0, elems)
    end
    if s == :joints1D
        return filter(elem -> elem.shape.class==JOINT1D_SHAPE, elems)
    end
    if s == :joints
        return filter(elem -> elem.shape.class==JOINT_SHAPE, elems)
    end
    if s == :nodes
        return get_nodes(elems)
    end
    if s == :ips
        return get_ips(elems)
    end
    error("Element getindex: Invalid symbol $s")
end

function getindex(elems::Array{Element,1}, cond::Expr)
    condm = fix_comparison_arrays(cond)
    funex = :( (x,y,z) -> false )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Element getindex: Invalid condition ", cond)
    end

    result = Array{Element}(0)
    for elem in elems
        coords = nodes_coords(elem.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        #if @static VERSION>v"0.6.0-rc1.0" ? Base.invokelatest(fun, x, y, z) : fun(x, y, z)
        if Base.invokelatest(fun, x, y, z)
            push!(result, elem) 
        end
    end
    return result
end

function getindex(elems::Array{Element,1}, cond::AbstractString) 
    # filter by tag
    if typeof(parse(cond)) == Symbol
        result = Array{Element}(0)
        for elem in elems
            if cond == elem.tag
                push!(result, elem) 
            end
        end
        return result
    end

    # filter by expression
    return getindex(elems, parse(cond))
end

