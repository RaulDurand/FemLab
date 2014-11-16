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
export Element
export set_mat, get_nodes, set_state, reset


abstract IpData

type Ip
    R::Array{Float64,1}
    w::Float64
    X::Array{Float64,1}
    id::Int
    #owner::Node
    data::IpData

    function Ip(R::Array, w::Float64)
        this = new(vec(R), w)
        this.X = Array(Float64, 3)
        this
    end
end


abstract Material


type Element
    shape ::ShapeType
    nodes ::Array{Node,1}
    ndim  ::Int
    tag   ::String
    id    ::Int
    ips   ::Array{Ip,1}
    active::Bool
    lnk_elems::Array{Element,1}
    mat::Material
    extra::Dict{Symbol,Any}

    function Element(shape, nodes, ndim, tag="")
        this     = new(shape, nodes, ndim, tag)
        this.ips = []
        this.extra = Dict{Symbol,Any}()
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

#L = [i+j; for i=1:10,j=1:10 ]
#L = [i+j  for i=1:10,j=1:10 if i>5 ]

function getindex(elems::Array{Element,1}, s::Symbol)
    if s == :solids
        cr = [ is_solid(elem.shape) for elem in elems]
    end
    if s == :lines
        cr = [ is_line(elem.shape) for elem in elems]
    end
    if s == :joints
        cr = [ is_joint(elem.shape) for elem in elems]
    end
    return elems[cr]
end

function getcoords(elem::Element)
    nnodes = length(elem.nodes)
    ndim   = elem.ndim
    return [ elem.nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end

function getconns(elem::Element)
    nnodes = length(elem.nodes)
    return [ node.id for node in elem.nodes ]
end

# dispatcher
config_dofs(elem::Element) = config_dofs(elem.mat, elem)

function set_mat(elem::Element, mm::Material, nips=0)
    ipc =  get_ip_coords(elem.shape, nips)
    nips = size(ipc,1)

    resize!(elem.ips, nips)
    elem.mat = mm
    for i=1:nips
        R = ipc[i,1:3]
        w = ipc[i,4]
        elem.ips[i] = Ip(R, w)
        elem.ips[i].id = i
        elem.ips[i].data = mm.new_ipdata(elem.ndim)
    end
    config_dofs(elem)
end

function set_mat(elems::Array{Element,1}, mm::Material)
    for elem in elems
        set_mat(elem, mm)
    end
end

function set_state(elems::Array{Element,1}; args...)
    for elem in elems
        for ip in elem.ips
            set_state(ip.data; args...)
        end
    end
end


function get_nodes(elems::Array{Element,1})
    nodes = Set{Node}()
    for elem in elems
        for node in elem.nodes
            push!(nodes, node)
        end
    end
    return [node for node in nodes]
end



