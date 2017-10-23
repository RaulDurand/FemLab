abstract type Monitor end

###############################################################################
# Monitor concrete structures


mutable struct NodeMonitor <: Monitor
    expr ::Expr
    table::DTable
    filename::String
    node::Node

    function NodeMonitor(node::Node, filename="")
        return new(:(), DTable(), filename, node)
    end

    #function NodeMonitor(nodes::Array{Node,1}, filename="")
        #if length(nodes)==0 
            #warn("NodeMonitor: Empty array of nodes to monitor.")
            #return new(:(), DTable(), filename)
        #end
        #length(nodes)>1 && warn("NodeMonitor: More than one node to monitor. Taking one.")
        #return new(:(), DTable(), filename, nodes[1])
    #end

    function NodeMonitor(expr::Expr, filename="")
        return new(expr, DTable(), filename)
    end
end


mutable struct NodesMonitor <: Monitor
    expr ::Expr
    book::DBook
    filename::String
    nodes::Array{Node,1}

    function NodesMonitor(nodes::Array{Node,1}, filename="")
        return new(:(), DBook(), filename, nodes)
    end

    function NodesMonitor(expr::Expr, filename="")
        return new(expr, DBook(), filename, [])
    end
end

abstract type FacetsMonitor <: Monitor end

mutable struct FacesMonitor <: FacetsMonitor
    expr ::Expr
    table::DTable
    filename::String
    nodes::Array{Node,1} # nodes from all selected facets

    function FacesMonitor(faces::Array{Face,1}, filename="")
        return new(:(), DTable(), filename, faces[:nodes])
    end

    function FacesMonitor(expr::Expr, filename="")
        return new(expr, DTable(), filename, [])
    end
end

mutable struct EdgesMonitor <: FacetsMonitor
    expr ::Expr
    table::DTable
    filename::String
    nodes::Array{Node,1} # nodes from all selected facets

    function EdgesMonitor(edges::Array{Edge,1}, filename="")
        return new(:(), DTable(), filename, edges[:nodes])
    end

    function EdgesMonitor(expr::Expr, filename="")
        return new(expr, DTable(), filename, [])
    end
end


mutable struct IpMonitor <: Monitor
    expr ::Expr
    table::DTable
    filename::String
    ip   ::Ip

    function IpMonitor(ip::Ip, filename="")
        return new(:(), DTable(), filename, ip)
    end

    #function IpMonitor(elem::Element, filename="")
        #if length(elem.ips)==0 
            #warn("IpMonitor: No ips found to monitor in array of elements. Check if material was defined for all elements")
            #return new(:(), DTable(), filename)
        #end
        #length(elem.ips)>1 && warn("IpMonitor: More than one ip to monitor with NodeMonitor. Taking one")
        #return new(:(), DTable(), filename, elem.ips[1])
    #end

    function IpMonitor(expr::Expr, filename="")
        return new(expr, DTable(), filename)
    end
end

mutable struct IpsMonitor <: Monitor
    expr ::Expr
    book ::DBook
    filename::String
    ips  ::Array{Ip,1}

    function IpsMonitor(ips::Array{Ip,1}, filename="")
        length(ips)==0 && warn("IpsMonitor: No ips found to monitor in array. Check if material was defined for all elements")
        return this = new(:(), DBook(), filename, ips)
    end

    #function IpsMonitor(elems::Array{Element,1}, filename="")
        #ips = get_ips(elems)
        #length(ips)==0 && warn("IpsMonitor: No ips found to monitor in array of elements. Check if material was defined for all elements")
        #return new(:(), DBook(), filename, ips)
    #end

    function IpsMonitor(expr::Expr, filename="")
        return new(expr, DBook(), filename, [])
    end
end

function save(monitor::Monitor, filename::AbstractString; verbose=true)
    if isdefined(monitor, :(table))
        save(monitor.table, filename, verbose=verbose)
    else
        save(monitor.book, filename, verbose=verbose)
    end
end

# For compatibility with older versions
const NodeTracker  = NodeMonitor
const NodesTracker = NodesMonitor
const IpTracker    = IpMonitor
const IpsTracker   = IpsMonitor
const FacesTracker = FacesMonitor
const EdgesTracker = EdgesMonitor


###############################################################################
# Functions to setup monitors


function setup_monitor!(monitor::NodeMonitor, domain)
    if monitor.expr != :() 
        nodes = domain.nodes[monitor.expr]
        if length(nodes) == 0
            warn("NodeMonitor: No nodes found for expression: $(monitor.expr)")
        elseif length(nodes) > 1
            warn("NodeMonitor: More than one node match expression: $(monitor.expr)")
        else
            monitor.node = nodes[1]
        end
    end
end


function setup_monitor!(monitor::NodesMonitor, domain)
    if monitor.expr != :() 
        monitor.nodes = domain.nodes[monitor.expr]
        length(monitor.nodes) == 0 && warn("NodesMonitor: No nodes found for expression: $(monitor.expr)")
    end
end


function setup_monitor!(monitor::FacesMonitor, domain)
    if monitor.expr != :() 
        edges = domain.faces[monitor.expr]
        length(edges) == 0 && warn("FacesMonitor: No faces found for expression: $(monitor.expr)")
        monitor.nodes = edges[:nodes]
    end
end

function setup_monitor!(monitor::EdgesMonitor, domain)
    if monitor.expr != :() 
        edges = domain.edges[monitor.expr]
        length(edges) == 0 && warn("EdgesMonitor: No edges found for expression: $(monitor.expr)")
        monitor.nodes = edges[:nodes]
    end
end


function setup_monitor!(monitor::IpMonitor, domain)
    if monitor.expr != :() 
        ips = domain.elems[:ips][monitor.expr]
        if length(ips) == 0
            warn("IpMonitor: No ips found for expression: $(monitor.expr)")
        elseif length(ips) > 1
            warn("IpMonitor: More than one ip match expression: $(monitor.expr)")
        else
            monitor.ip = ips[1]
        end
    end
end


function setup_monitor!(monitor::IpsMonitor, domain)
    if monitor.expr != :() 
        monitor.ips = domain.elems[:ips][monitor.expr]
        length(monitor.ips) == 0 && warn("IpsMonitor: No ips found for expression: $(monitor.expr)")
    end
end


###############################################################################
# Functions to update monitors


function update_monitor_file(monitor::Monitor)
    if monitor.filename != ""
        try   save(monitor, monitor.filename, verbose=false)
        catch warn("Problem writing file ", monitor.filename)
        end
    end
end


function update_monitor!(monitor::NodeMonitor)
    !isdefined(monitor, :node) && return

    vals = node_vals(monitor.node)
    push!(monitor.table, vals)

    update_monitor_file(monitor)
end

function update_monitor!(monitor::NodesMonitor)
    table = DTable()
    for node in monitor.nodes
        vals = node_vals(node)
        push!(table, vals)
    end
    push!(monitor.book, table)

    update_monitor_file(monitor)
end

function update_monitor!(monitor::FacetsMonitor)
    # update data
    tableU = DTable()
    tableF = DTable()
    for node in monitor.nodes
        valsU  = Dict( dof.sU::Symbol => dof.U::Float64 for dof in node.dofs )
        valsF  = Dict( dof.sF::Symbol => dof.F::Float64 for dof in node.dofs )
        push!(tableF, valsF)
        push!(tableU, valsU)
    end

    valsU = Dict( key => mean(tableU[key]) for key in keys(tableU) ) # gets the average of essential values
    valsF = Dict( key => sum(tableF[key])  for key in keys(tableF) ) # gets the sum for each component
    vals  = merge(valsU, valsF)
    push!(monitor.table, vals)

    update_monitor_file(monitor)
end

function update_monitor!(monitor::IpMonitor)
    !isdefined(monitor, :ip) && return

    vals = ip_state_vals(monitor.ip.owner.mat, monitor.ip.data)
    push!(monitor.table, vals)

    update_monitor_file(monitor)
end

function update_monitor!(monitor::IpsMonitor)
    table = DTable()
    for ip in monitor.ips
        vals = ip_vals(ip)  # includes ip global coordinates
        push!(table, vals)
    end
    push!(monitor.book, table)

    update_monitor_file(monitor)
end
