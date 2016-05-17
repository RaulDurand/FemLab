abstract Tracker

export NodeTracker, NodesTracker, IpTracker, IpsTracker

type NodeTracker <: Tracker
    table::DTable
    node::Node

    function NodeTracker(node::Node)
        return new(DTable(), node )
    end

    function NodeTracker(nodes::Array{Node,1})
        return new(DTable(), nodes[1])
    end

    function NodeTracker(domain, expr::Expr) # TODO: verify if this syntax is convenient
        nodes = domain.nodes[expr]
        @assert length(nodes)==1
        return new(DTable(), nodes[1])
    end
end

type NodesTracker <: Tracker
    book::DBook
    nodes::Array{Node,1}

    function NodesTracker(nodes::Array{Node,1})
        return new(DBook(), nodes[:])
    end

    function NodesTracker(domain, expr::Expr)
        nodes = domain.nodes[expr]
        @assert length(nodes)>0
        return new(DBook(), nodes)
    end
end

type IpTracker <: Tracker
    table::DTable
    ip   ::Ip

    function IpTracker(ip::Ip)
        return new(DTable(), ip)
    end

    function IpTracker(elem::Element)
        @assert length(elem.ips)>0
        return new(DTable(), elem.ips[1])
    end

    function IpTracker(domain, expr::Expr)
        ips = get_ips(domain.elems)[expr]
        @assert length(ips)==1
        return new(DTable(), ips[1])
    end
end

type IpsTracker <: Tracker
    book::DBook
    ips  ::Array{Ip,1}

    function IpsTracker(ips::Array{Ip,1})
        return this = new(DBook(), ips)
    end

    function IpsTracker(elems::Array{Element,1})
        ips = get_ips(elems)
        @assert length(ips)>0
        return new(DBook(), ips)
    end

    function IpsTracker(domain, expr::Expr)
        ips = get_ips(domain.elems)[expr]
        @assert length(ips)>0
        return new(DBook(), ips)
    end
end

function save(trk::Tracker, filename::AbstractString; verbose=true, format="dat")
    if typeof(trk) in (NodeTracker, IpTracker)
        save(trk.table, filename, verbose=verbose, format=format)
    else
        save(trk.book , filename, verbose=verbose, format=format)
    end
end
