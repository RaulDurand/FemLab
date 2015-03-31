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
export set_mat, get_nodes, get_ips, set_state, reset
export getvals
export read_prms


# Abstract type for IP data
# =========================

abstract IpData


#  Ip
# ====

type Ip
    R    ::Array{Float64,1}
    w    ::Float64
    X    ::Array{Float64,1}
    id   ::Int
    owner::Any    # Element
    data ::IpData
    data0::IpData # Backup

    function Ip(R::Array, w::Float64)
        this   = new(vec(R), w)
        this.X = Array(Float64, 3)
        this.owner = nothing
        this
    end
end


# Get ip values in a dictionary
function getvals(ip::Ip)
    coords = [ :x => ip.X[1], :y => ip.X[1], :z => ip.X[3] ]
    vals   = getvals(ip.owner.mat, ip.data)
    return merge(coords, vals)
end

# Index operator for a ip collection using expression
function getindex(ips::Array{Ip,1}, cond::Expr) 
    condm = fix_comparison_scalar(cond)
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Ip getindex: Invalid condition ", cond)
    end

    result = Array(Ip, 0)
    for ip in ips
        if fun(ip.X[1], ip.X[2], ip.X[3])
            push!(result, ip)
        end
    end
    return result
end

getindex(ips::Array{Ip,1}, cond::String) = getindex(ips, parse(cond))



# Abstract tyep for material
# ==========================

abstract Material



# Type Element
# ============

type Element
    shape ::ShapeType
    nodes ::Array{Node,1}
    ndim  ::Int
    tag   ::String
    id    ::Int
    ips   ::Array{Ip,1}
    active::Bool
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

# Index operator for a collection of elements
function getindex(elems::Array{Element,1}, s::Symbol)
    if s == :solids
        cr = [ is_solid(elem.shape) for elem in elems]
        return elems[cr]
    end
    if s == :lines
        cr = [ is_line(elem.shape) for elem in elems]
        return elems[cr]
    end
    if s == :joints
        cr = [ is_joint(elem.shape) for elem in elems]
        return elems[cr]
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
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Element getindex: Invalid condition ", cond)
    end

    result = Array(Element,0)
    for elem in elems
        coords = getcoords(elem.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        if fun(x, y, z)
            push!(result, elem) 
        end
    end
    return result
end

function getindex(elems::Array{Element,1}, cond::String) 
    if typeof(parse(cond)) == Symbol
        result = Array(Element,0)
        for elem in elems
            if cond == elem.tag
                push!(result, elem) 
            end
        end
        return result
    end

    return getindex(elems, parse(cond))
end

# Get the element coordinates matrix
function getcoords(elem::Element)
    nnodes = length(elem.nodes)
    ndim   = elem.ndim
    return [ elem.nodes[i].X[j] for i=1:nnodes, j=1:ndim]
end

function getconns(elem::Element)
    nnodes = length(elem.nodes)
    return [ node.id for node in elem.nodes ]
end


# Dispatcher: Configure all dofs in an element according to material
config_dofs(elem::Element) = config_dofs(elem.mat, elem)


# Define material properties for an element
function set_mat(elem::Element, mm::Material; nips::Int64=0)
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
        elem.ips[i].owner = elem
    end

    # finding ips global coordinates
    C     = getcoords(elem)
    shape = elem.shape
    if is_joint(shape)
        C     = getcoords(elem.extra[:bar])
        shape = elem.extra[:bar].shape
    end

    for ip in elem.ips
        N = shape_func(shape, ip.R)
        ip.X = C'*N
    end

    # configure degrees of freedom
    config_dofs(elem)
end


# Define material properties for a collection of elements
function set_mat(elems::Array{Element,1}, mm::Material; nips::Int64=0)
    for elem in elems
        set_mat(elem, mm, nips=nips)
    end
end


# Reads material parameters from a json file
function read_prms(filename::String)

    # read file
    file = open(filename, "r")
    data = JSON.parse(file)
    close(file)

    # parse materials
    mats_prms = Dict{String, Any}()
    for d in data
        name = d["name"]
        keys = d["prms"]
        vals = d["vals"]
        prms = (Symbol => Float64)[ symbol(k) => v for (k,v) in zip(keys, vals) ]
        mats_prms[name] = prms
    end

    return mats_prms
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

