abstract Monitor

import Base.getindex
export NodeTracker, NodesTracker, IpTracker, IpsTracker, FacesTracker

type NodeTracker <: Monitor
    table::DTable
    node::Node
    filename::String

    function NodeTracker(node::Node, filename="")
        return new(DTable(), node, filename )
    end

    function NodeTracker(nodes::Array{Node,1}, filename="")
        return new(DTable(), nodes[1], filename)
    end

    function NodeTracker(domain, expr::Expr, filename="") # TODO: verify if this syntax is convenient
        nodes = domain.nodes[expr]
        @assert length(nodes)==1
        return new(DTable(), nodes[1], filename)
    end
end

type NodesTracker <: Monitor
    book::DBook
    nodes::Array{Node,1}
    filename::String

    function NodesTracker(nodes::Array{Node,1}, filename="")
        return new(DBook(), nodes[:], filename)
    end

    function NodesTracker(domain, expr::Expr, filename="")
        nodes = domain.nodes[expr]
        @assert length(nodes)>0
        return new(DBook(), nodes, filename)
    end
end

type FacesTracker <: Monitor
    table::DTable
    nodes::Array{Node,1} # nodes from all selected faces
    filename::String

    function FacesTracker(faces::Array{Face,1}, filename="")
        return new(DTable(), faces[:nodes], filename)
    end

    function FacesTracker(domain, expr::Expr, filename="")
        faces = domain.faces[expr]
        @assert length(faces)>0
        return new(DTable(), faces[:nodes], filename)
    end
end

type IpTracker <: Monitor
    table::DTable
    ip   ::Ip
    filename::String

    function IpTracker(ip::Ip, filename="")
        return new(DTable(), ip, filename)
    end

    function IpTracker(elem::Element, filename="")
        @assert length(elem.ips)>0
        return new(DTable(), elem.ips[1], filename)
    end

    function IpTracker(domain, expr::Expr, filename="")
        ips = get_ips(domain.elems)[expr]
        @assert length(ips)==1
        return new(DTable(), ips[1], filename)
    end
end

type IpsTracker <: Monitor
    book::DBook
    ips  ::Array{Ip,1}
    filename::String

    function IpsTracker(ips::Array{Ip,1}, filename="")
        return this = new(DBook(), ips, filename)
    end

    function IpsTracker(elems::Array{Element,1}, filename="")
        ips = get_ips(elems)
        @assert length(ips)>0
        return new(DBook(), ips, filename)
    end

    function IpsTracker(domain, expr::Expr, filename="")
        ips = get_ips(domain.elems)[expr]
        @assert length(ips)>0
        return new(DBook(), ips, filename)
    end
end

function save(monit::Monitor, filename::AbstractString; verbose=true, format="dat")
    if typeof(monit) in (NodeTracker, IpTracker, FacesTracker)
        save(monit.table, filename, verbose=verbose, format=format)
    else
        save(monit.book , filename, verbose=verbose, format=format)
    end
end
