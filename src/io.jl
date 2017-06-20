import Base.show


macro show_array_function(datatype)
    return quote
        function $(esc(:show))(io::IO, array::Array{$(datatype),1})
            n = length(array)
            maxn = 12
            half = div(maxn,2)
            idx = n<=maxn ? [1:n;] : [1:half; n-half+1:n]
            for i in idx
                println(io, array[i])
                if n>maxn && i==half
                    println(io, "⋮")
                end
            end
            return nothing
        end
    end
end



# Show basic information of a dof using print
function show(io::IO, dof::Dof)
    @printf io "Dof( eq_id=%s  %s=%-8.4f    %s=%-8.4f    bry_%s=%-8.4f   bry_%s=%-8.4f  presc=%s)" dof.eq_id   dof.sU dof.U dof.sF dof.F   dof.sU dof.bryU dof.sF dof.bryF dof.prescU
end

@show_array_function Dof



# Show basic information of a node using print
function show(io::IO, node::Node)
    X = node.X
    sdofs = "Dof"*string( [dof.sU for dof in node.dofs] )
    @printf io "Node( id=%s X=[%-12.4f, %-12.4f, %-12.4f] tag=%s, dofs=%s )" node.id X[1] X[2] X[3] node.tag sdofs
end

@show_array_function Node



# Show basic information of an ip using print
function show(io::IO, ip::Ip)
    datas = isdefined(ip, :data) ? "$(typeof(ip.data))()" : ""
    @printf io "Ip( R=%s    X=%s    data=%s )" rpad(ip.R, 25, " ") rpad(ip.X, 25, " ") datas
end

@show_array_function Ip



# Show basic information of an element using print
function show(io::IO, elem::Element)
    snodes = "Node"*string([n.id for n in elem.nodes])
    sips   = length(elem.ips)==0 ? "Ip[]" : "Ip[...]"
    @printf io "Element( id=%s  shape=%s  tag=%s  active=%s  nodes=%s  ips=%s  mat=%s )" elem.id elem.shape.name elem.tag elem.active snodes sips typeof(elem.mat)
end

@show_array_function Element



function show(io::IO, facet::Facet)
    snodes = "Node"*string([n.id for n in facet.nodes])
    soelem = isdefined(facet, :oelem) ? "$(facet.oelem.id)" : ""
    @printf io "Facet( shape=%s  nodes=%s  oelem=%s  isedge=%s )" facet.shape.name snodes soelem facet.isedge
end

@show_array_function Facet


function show(io::IO, bc::NodesBC)
    nnodes = length(bc.nodes)
    if nnodes>0
        println(io, "NodeBC( nodes=Node[...] total $nnodes   conds=($(bc.conds))")
    else
        println(io, "NodeBC( expr=($(bc.expr))   conds=$(bc.conds)")
    end
end

function show(io::IO, bc::FacesBC)
    nfaces = length(bc.faces)
    if nfaces>0
        println(io, "FaceBC( faces=Face[...] total $nfaces   conds=($(bc.conds))")
    else
        println(io, "FaceBC( expr=($(bc.expr))   conds=$(bc.conds)")
    end
end

function show(io::IO, bc::EdgesBC)
    nedges = length(bc.edges)
    if nedges>0
        println(io, "FaceBC( edges=Edge[...] total $nedges   conds=($(bc.conds))")
    else
        println(io, "FaceBC( expr=($(bc.expr))   conds=$(bc.conds)")
    end
end

@show_array_function BC


function show(io::IO, dom::Domain)
    nnodes = length(dom.nodes)
    nelems = length(dom.elems)
    out = """Domain(
    ndim = $(dom.ndim)
    nodes = $(nnodes>0 ? "[...]" : "[]")   total $nnodes
    elems = $(nelems>0 ? "[...]" : "[]")   total $nelems
    )"""
    println(io, out)
end

