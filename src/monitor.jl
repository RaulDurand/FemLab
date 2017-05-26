abstract Monitor

type NodeMonitor <: Monitor
    table::DTable
    node::Node
    filename::String

    function NodeMonitor(node::Node, filename="")
        return new(DTable(), node, filename )
    end

    function NodeMonitor(nodes::Array{Node,1}, filename="")
        return new(DTable(), nodes[1], filename)
    end

    function NodeMonitor(domain, expr::Expr, filename="") # TODO: verify if this syntax is convenient
        nodes = domain.nodes[expr]
        @assert length(nodes)==1
        return new(DTable(), nodes[1], filename)
    end
end

type NodesMonitor <: Monitor
    book::DBook
    nodes::Array{Node,1}
    filename::String

    function NodesMonitor(nodes::Array{Node,1}, filename="")
        return new(DBook(), nodes[:], filename)
    end

    function NodesMonitor(domain, expr::Expr, filename="")
        nodes = domain.nodes[expr]
        @assert length(nodes)>0
        return new(DBook(), nodes, filename)
    end
end

type FacesMonitor <: Monitor # TODO: change to FacetsMonitor
    table::DTable
    nodes::Array{Node,1} # nodes from all selected faces
    filename::String

    function FacesMonitor(faces::Array{Face,1}, filename="")
        return new(DTable(), faces[:nodes], filename)
    end

    function FacesMonitor(domain, expr::Expr, filename="")
        faces = domain.faces[expr]
        @assert length(faces)>0
        return new(DTable(), faces[:nodes], filename)
    end
end

EdgesMonitor = FacesMonitor

type IpMonitor <: Monitor
    table::DTable
    ip   ::Ip
    filename::String

    function IpMonitor(ip::Ip, filename="")
        return new(DTable(), ip, filename)
    end

    function IpMonitor(elem::Element, filename="")
        @assert length(elem.ips)>0
        return new(DTable(), elem.ips[1], filename)
    end

    function IpMonitor(domain, expr::Expr, filename="")
        ips = get_ips(domain.elems)[expr]
        @assert length(ips)==1
        return new(DTable(), ips[1], filename)
    end
end

type IpsMonitor <: Monitor
    book::DBook
    ips  ::Array{Ip,1}
    filename::String

    function IpsMonitor(ips::Array{Ip,1}, filename="")
        return this = new(DBook(), ips, filename)
    end

    function IpsMonitor(elems::Array{Element,1}, filename="")
        ips = get_ips(elems)
        @assert length(ips)>0
        return new(DBook(), ips, filename)
    end

    function IpsMonitor(domain, expr::Expr, filename="")
        ips = get_ips(domain.elems)[expr]
        @assert length(ips)>0
        return new(DBook(), ips, filename)
    end
end

function save(monit::Monitor, filename::AbstractString; verbose=true, format="dat")
    if typeof(monit) in (NodeMonitor, IpMonitor, FacesMonitor)
        save(monit.table, filename, verbose=verbose, format=format)
    else
        save(monit.book, filename, verbose=verbose, format=format)
    end
end

const NodeTracker  = NodeMonitor
const NodesTracker = NodesMonitor
const IpTracker    = IpMonitor
const IpsTracker   = IpsMonitor
const FacesTracker = FacesMonitor
const FacesTracker = EdgesMonitor
