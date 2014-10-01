
export set_bc, sort, reset


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
    data_table::Array
    function Node(X::Array{Float64,1}; tag="", id=-1)
        this = new(X, tag, id)
        this.dofs = []
        this.dofdict = Dict()
        this.data_table = []
        this
    end
    function Node(point::Point; tag="", id=-1)
        Node([point.x, point.y, point.z], tag=tag, id=id)
    end
end

import Base.reset
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

function get_vals(node::Node)
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
function substitute!(expr::Expr, syms::Tuple, vals::Tuple)
    for (i,s) in enumerate(expr.args)
        if typeof(s)==Expr; substitute!(s, syms, vals); continue end
        for (sym,val) in zip(syms,vals)
            if s==sym; expr.args[i] = val end
        end
    end
end

function getindex(nodes::Array{Node,1}, cond::Expr) 
    result = Array(Node, 0)
    for node in nodes
        mcond = copy(cond)
        substitute!(mcond, (:x, :y, :z), (node.X[1], node.X[2], node.X[3]) )
        if eval(mcond); push!(result, node) end
    end
    return result
end

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

import Base.max
function max(nodes::Array{Node,1}, d::Symbol) 
    idx = findfisrt((:x, :y, :z), d)
    max([node.x[idx] for node in nodes]...)
end

import Base.min
function min(nodes::Array{Node,1}, d::Symbol) 
    idx = findfisrt((:x, :y, :z), d)
    min([node.x[idx] for node in nodes]...)
end

import Base.sort
function sort(nodes::Array{Node,1}, d::Symbol) 
    idx = findfirst((:x, :y, :z), d)
    idxs = sortperm([node.X[idx] for node in nodes])
    return nodes[idxs]
end

