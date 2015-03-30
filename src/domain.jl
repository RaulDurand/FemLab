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

export Domain
export track
export get_node

type Face
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union(Element,Nothing)
    isedge::Bool
    function Face(shape, nodes, ndim)
        this = new(shape, nodes, ndim)
        this.oelem  = nothing
        this.isedge = false
        return this
    end
end

typealias Edge Face


function fix_comparison_arrays(expr::Expr)
    mexpr = copy(expr) # expression to be modified
    tol = 1e-6

    fix_comp = function(expr::Expr)
        symb = expr.args[2]
        a = expr.args[1]
        b = expr.args[3]
        if symb == :(==)
            return :(maximum(abs($a-$b)) < $tol)
        end
        if symb == :(>=)
            return :(minimum($a) > maximum($b) - $tol)
        end
        if symb == :(>)
            return :(minimum($a) > maximum($b) + $tol)
        end
        if symb == :(<=)
            return :(maximum($a) < minimum($b) + $tol)
        end
        if symb == :(<)
            return :(maximum($a) < minimum($b) - $tol)
        end
        return expr
    end

    if mexpr.head == :comparison
        return fix_comp(mexpr)
    else
        for (i,arg) in enumerate(mexpr.args)
            if typeof(arg)!=Expr; continue end
            if arg.head == :comparison
                mexpr.args[i] = fix_comp(arg)
            else
                mexpr.args[i] = fix_comparison_arrays(arg)
            end
        end
    end
    return mexpr
end

function getindex(faces::Array{Face,1}, cond::Expr)
    #condm = subs_equal_by_approx(cond)
    condm = fix_comparison_arrays(cond)
    #@show condm
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Faces getindex: Invalid condition ", cond)
    end

    result = Array(Face,0)
    for face in faces
        coords = getcoords(face.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        if fun(x, y, z)
            push!(result, face) 
        end
    end
    return result
end

getindex(faces::Array{Face,1}, cond::String) = getindex(faces, parse(cond))

function set_bc(face::Face; args...)
    oelem = face.oelem # owner element
    if oelem==nothing; error("Face with no owner element") end

    for (key,val) in args
        set_facet_bc(oelem.mat, oelem, face, key, float(val))
    end
end

function set_bc(faces::Array{Face,1}; args...)
    if length(faces)==0; println(RED, "Warning, applying boundary conditions to empty array of faces", DEFAULT) end
    for face in faces
        set_bc(face; args...)
    end
end


type Domain
    ndim ::Int
    nodes::Array{Node,1}
    elems::Array{Element,1}
    faces::Array{Face,1}
    edges::Array{Edge,1}
    filekey::String

    trk_nodes     ::Array{(Node, DTable), 1}
    trk_ips       ::Array{(Ip  , DTable), 1}
    trk_list_nodes::Array{(Array{Node,1}, DBook), 1}
    trk_list_ips  ::Array{(Array{Ip  ,1}, DBook), 1}
    function Domain(mesh::Mesh; filekey="out")
        this = new()
        this.filekey = filekey

        load_mesh(this, mesh)
        this.trk_nodes = []
        this.trk_ips   = []
        this.trk_list_nodes = []
        this.trk_list_ips   = []
        return this
    end
end


function load_mesh(dom::Domain, mesh::Mesh)
    ndim = dom.ndim = mesh.ndim

    # Setting nodes
    dom.nodes = [ Node(point, id=i) for (i,point) in enumerate(mesh.points)]

    # Setting elements
    dom.elems = Array(Element,0)
    for (i,cell) in enumerate(mesh.cells)
        conn = [ p.id for p in cell.points ]
        elem = Element(cell.shape, dom.nodes[conn], ndim)
        elem.id = i
        push!(dom.elems, elem)
    end

    # Setting faces
    dom.faces = Array(Face,0)
    for cell in mesh.faces
        conn = [ p.id for p in cell.points ]
        face = Face(cell.shape, dom.nodes[conn], ndim)
        face.oelem = dom.elems[cell.ocell.id]
        push!(dom.faces, face)
    end

    # Setting edges
    dom.edges = Array(Edge,0)
    for cell in mesh.edges
        conn = [ p.id for p in cell.points ]
        edge = Face(cell.shape, dom.nodes[conn], ndim)
        edge.oelem = dom.elems[cell.ocell.id]
        push!(dom.edges, edge)
    end

    # Setting embeddeds
    edict = Dict{Uint64, Element}()
    for elem in dom.elems
        #hs = hash(getcoords(elem))
        hs = hash(getconns(elem))
        edict[hs] = elem
    end
    for elem in dom.elems
        if elem.shape in (LINK2, LINK3)
            nbnodes = elem.shape==LINK2? 2 : 3
            conns = getconns(elem)
            hs_hook  = hash(conns[1:end-nbnodes])
            hs_truss = hash(conns[end-nbnodes+1:end])
            elem.extra[:hook] = edict[hs_hook]
            elem.extra[:bar ] = edict[hs_truss]
        end
    end
    
end

#macro 


function calc_nodal_vals(dom::Domain)
    # Get incidence matrix (shares) (fast)
    np = length(points)
    Shares = [ Element[] for i=1:np]
    for elem in dom.elems
        for node in elem.nodes
            push!(Shares[node.id], elem)
        end
    end

    # Interpolating solid elements
    #for node in dom.nodes
        #patch = @list(elem, elem in node.shares, if issolid(elem) end)
        #patch = node.shares[ map(is_solid, node.shares) ]
    #end
end

function make_backup(elems::Array{Element,1})
    for elem in elems
        for ip in elem.ips
            ip.data0 = deepcopy(ip.data)
        end
    end
end

function node_and_elem_vals(nodes::Array{Node,1}, elems::Array{Element,1})

    # Calculate max number of node and elem labes
    matset  = Set{Material}()
    nlabels = Set{Symbol}()
    elabels = Set{Symbol}()
    for elem in elems
        if !(elem.mat in matset)
            # getting values from current element
            node_vals, elem_vals = node_and_elem_vals(elem.mat, elem)
            union!(nlabels, keys(node_vals))
            union!(elabels, keys(elem_vals))
            push!(matset, elem.mat)
        end
    end

    nlabels_idx = [ key=>i for (i,key) in enumerate(nlabels) ]
    elabels_idx = [ key=>i for (i,key) in enumerate(elabels) ]
    nncomps     = length(nlabels)
    necomps     = length(elabels)

    nnodes     = length(nodes)
    nelems     = length(elems)
    Node_vals  = zeros(Float32, nnodes, nncomps) # float32 is important for paraview
    Node_reps  = zeros(Int64  , nnodes, nncomps)
    Elem_vals  = zeros(Float32, nelems, necomps) # float32 is important for paraview

    # In elements
    for elem in elems
        # getting values from current element
        node_vals, elem_vals = node_and_elem_vals(elem.mat, elem)

        # filling nodal values and node repetitions
        for (key, values) in node_vals
            idx = nlabels_idx[key]
            for (i, node) in enumerate(elem.nodes)
                Node_reps[node.id, idx] += 1
                Node_vals[node.id, idx] += values[i]
            end
        end

        # filling element values
        for (key,value) in elem_vals
            idx = elabels_idx[key]
            Elem_vals[elem.id, idx] = value
        end
    end

    # averaging nodal values
    for i = 1:nnodes
        for j=1:length(nlabels_idx)
            if Node_reps[i,j]>0
                Node_vals[i,j] /= Node_reps[i,j]
            end
        end
    end

    # Filling nodal labels list
    node_keys = [ nlabels... ]
    elem_keys = [ elabels... ]

    return Node_vals, node_keys, Elem_vals, elem_keys
end


function node_and_elem_vals2(nodes::Array{Node,1}, elems::Array{Element,1})

    # Calculate max number of node and elem labes
    matset  = Set{Material}()
    nlabels = Set{Symbol}()
    elabels = Set{Symbol}()
    for elem in elems
        if !(elem.mat in matset)
            # getting values from current element
            node_vals, elem_vals = node_and_elem_vals(elem.mat, elem)
            union!(nlabels, keys(node_vals))
            union!(elabels, keys(elem_vals))
            push!(matset, elem.mat)
        end
    end

    nlabels_idx = [ key=>i for (i,key) in enumerate(nlabels) ]
    elabels_idx = [ key=>i for (i,key) in enumerate(elabels) ]
    nncomps     = length(nlabels)
    necomps     = length(elabels)

    nnodes     = length(nodes)
    nelems     = length(elems)
    Node_vals  = zeros(Float32, nnodes, nncomps) # float32 is important for paraview
    Node_reps  = zeros(Int64  , nnodes, nncomps)

    # In elements
    for elem in elems
        # getting values from current element
        node_vals, elem_vals = node_and_elem_vals(elem.mat, elem)

        # filling nodal values and node repetitions
        for (key, values) in node_vals
            idx = nlabels_idx[key]
            for (i, node) in enumerate(elem.nodes)
                Node_reps[node.id, idx] += 1
                Node_vals[node.id, idx] += values[i]
            end
        end

        # filling element values
        #for (key,value) in elem_vals
            #idx = elabels_idx[key]
            #Elem_vals[elem.id, idx] = value
        #end
    end

    # averaging nodal values
    for i = 1:nnodes
        for j=1:length(nlabels_idx)
            if Node_reps[i,j]>0
                Node_vals[i,j] /= Node_reps[i,j]
            end
        end
    end

    # Get element vals
    Elem_vals  = DTable()

    # Filling nodal labels list
    node_keys = [ nlabels... ]
    elem_keys = [ elabels... ]

    return Node_vals, node_keys, Elem_vals, elem_keys
end

function get_node(dom::Domain, coord::Array{Float64,1})
    X = [ coord, 0.0 ][1:3]
    tol     = 1.0e-8
    for node in dom.nodes
        if norm(X-node.X) < tol
            return node
        end
    end
    return nothing
end

function track(dom::Domain, node::Node)
    new_table = DTable()
    push!(dom.trk_nodes, (node, new_table))
    return new_table
end

function track(dom::Domain, ip::Ip)
    new_table = DTable()
    push!(dom.trk_ips, (ip, new_table))
    return new_table
end

function track(dom::Domain, elem::Element)
    new_table = DTable()
    push!(dom.trk_ips, (elem.ips[1], new_table))
    return new_table
end

function track(dom::Domain, node_list::Array{Node,1})
    new_book = DBook()
    push!(dom.trk_list_nodes, (node_list, new_book))
    return new_book
end

function track(dom::Domain, ip_list::Array{Ip,1})
    new_book = DBook()
    push!(dom.trk_list_ips, (ip_list, new_book))
    return new_book
end

function tracking(dom::Domain)
    # tracking nodes
    for (node, table) in dom.trk_nodes
        vals = getvals(node)
        push!(table, vals)
    end

    # tracking ips
    for (ip, table) in dom.trk_ips
        vals = getvals(ip.owner.mat, ip.data)
        push!(table, vals)
    end

    # tracking list of nodes
    for (nodes, book) in dom.trk_list_nodes
        table = DTable()
        for node in nodes
            vals = getvals(node)
            push!(table, vals)
        end
        push!(book, table)
    end

    # tracking list of ips
    for (ips, book) in dom.trk_list_ips
        table = DTable()

        for ip in ips
            vals = getvals(ip)
            push!(table, vals)
        end
        push!(book, table)
    end
end



function save(dom::Domain, filename::String; verbose=true, save_ips=false)
    # Saves the dom information in vtk format
    nnodes = length(dom.nodes)
    nelems  = length(dom.elems)

    # Number of total connectivities
    nconns = 0
    for elem in dom.elems
        nconns += 1 + length(elem.nodes)
    end

    # Open filename
    f = open(filename, "w")

    println(f, "# vtk DataFile Version 3.0")
    println(f, "pyfem output ")
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", nnodes, " float64")

    # Write nodes
    for node in dom.nodes
        @printf f "%23.15e %23.15e %23.15e \n" node.X[1] node.X[2] node.X[3]
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ",nelems, " ", nconns)
    for elem in dom.elems
        print(f, length(elem.nodes), " ")
        for node in elem.nodes
            print(f, node.id-1, " ")
        end
        println(f)
    end
    println(f)

    # Write elem types
    println(f, "CELL_TYPES ", nelems)
    for elem in dom.elems
        println(f, get_vtk_type(elem.shape))
    end
    println(f)

    # check if all elements have material defined
    has_data = all([ isdefined(elem, :mat) for elem in dom.elems])

    if has_data
        # Get node and elem values
        node_vals, node_labels, elem_vals, elem_labels = node_and_elem_vals(dom.nodes, dom.elems)
        nncomps = length(node_labels)
        necomps = length(elem_labels)
    else
        close(f)
        return
    end

    # Write point data
    println(f, "POINT_DATA ", nnodes)

    # Write vectors
    if :ux in keys(dom.nodes[1].dofdict)
        println(f, "VECTORS ", "|u| float64")
        for node in dom.nodes
            @printf f "%23.10e"   node.dofdict[:ux].U
            @printf f "%23.10e"   node.dofdict[:uy].U
            @printf f "%23.10e\n" dom.ndim==3?node.dofdict[:uz].U:0.0
        end
    end
    println(f, )

    # Write nodal scalar data
    for i=1:nncomps
        println(f, "SCALARS ", node_labels[i], " float64 1")
        println(f, "LOOKUP_TABLE default")
        for j in 1:nnodes
            @printf f "%23.10e"  node_vals[j,i]
        end
        println(f, )
    end


    # Write element data
    println(f, "CELL_DATA ", nelems)
    for i=1:necomps
        println(f, "SCALARS ", elem_labels[i], " float64 1")
        println(f, "LOOKUP_TABLE default")
        for j in 1:nelems  #naelems
            e_idx = dom.elems[j].id  #aelems
            @printf f "%23.10e"  elem_vals[e_idx,i]
        end
        println(f, )
    end

    # Write elem type
    println(f, "SCALARS cell_type int 1")
    println(f, "LOOKUP_TABLE default")
    for elem in dom.elems
        println(f, elem.shape, " ")
    end
    println(f, )

    close(f)

    if verbose
        println(GREEN, "  file $filename written (Domain)", DEFAULT)
    end

    # save ip information as vtk
    if has_data && save_ips
        basename, ext = splitext(filename)
        save_dom_ips(dom, basename*"-ip"*ext, verbose)
    end

end


function save_dom_ips(dom::Domain, filename::String, verbose=true)
    # Saves ips information from domain as a vtk file

    # Get all ips
    ips = Ip[]
    for elem in dom.elems
        for ip in elem.ips
            push!(ips, ip)
        end
    end

    nips   = length(ips)
    ncells = nips

    # Number of total connectivities
    nconns = nips*2 

    # Open filename
    f = open(filename, "w")

    println(f, "# vtk DataFile Version 3.0")
    println(f, "FemLab output ")
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", nips, " float64")

    # Write ip points
    for ip in ips
        x = ip.X[1]
        y = ip.X[2]
        z = length(ip.X)>2? ip.X[3] : 0.0
        @printf f "%23.15e %23.15e %23.15e \n" x y z
    end
    println(f)

    # Write ip connectivities (vertex)
    println(f, "CELLS ", ncells, " ", nconns)
    for (i,ip) in enumerate(ips)
        println(f, "1 $(i-1)", )
    end
    println(f)

    # Ips cell types
    println(f, "CELL_TYPES ", nips)
    vtk_vertex = 1
    for ip in ips
        print(f, vtk_vertex, " ")
    end
    println(f)

    # Get values at ips
    table = DTable()
    for ip in ips
        push!(table, getvals(ip))
    end

    # Write point data
    println(f, "POINT_DATA ", nips)

    # Write ip scalar data TODO
    for field in table.header
        println(f, "SCALARS ", field, " float64 1")
        println(f, "LOOKUP_TABLE default")
        field_data = table.dict[field]
        for i in 1:nips
            @printf f "%23.10e" field_data[i]
        end
        println(f, )
    end

    close(f)

    if verbose
        println(GREEN, "  file $filename written (Domain)", DEFAULT)
    end
end


function save2(dom::Domain, filename::String; verbose=true)
    # Saves the dom information in vtk format
    nnodes = length(dom.nodes)
    nelems = length(dom.elems)
    npoints = nnodes

    # check if all elements have material defined
    has_data = all([ isdefined(elem, :mat) for elem in dom.elems])

    ips = Ip[]
    if has_data
        # Get all ips
        for elem in dom.elems
            for ip in elem.ips
                push!(ips, ip)
            end
        end
    end
    nips = length(ips)

    # Number of total connectivities
    nconns = 0
    for elem in dom.elems
        nconns += 1 + length(elem.nodes)
    end

    # Add integration points
    nconns  += nips*2
    npoints += nips
    ncells   = nelems + nips

    # Open filename
    f = open(filename, "w")

    println(f, "# vtk DataFile Version 3.0")
    println(f, "FemLab output ")
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", npoints, " float64")

    # Write nodes
    for node in dom.nodes
        @printf f "%23.15e %23.15e %23.15e \n" node.X[1] node.X[2] node.X[3]
    end
    # Write ip points
    for ip in ips
        @printf f "%23.15e %23.15e %23.15e \n" ip.X[1] ip.X[2] ip.X[3]
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ", ncells, " ", nconns)
    for elem in dom.elems
        print(f, length(elem.nodes), " ")
        for node in elem.nodes
            print(f, node.id-1, " ")
        end
        println(f)
    end
    # Write ip connectivities (vertex)
    for (i,ip) in enumerate(ips)
        println(f, "1 $(i-1)", )
    end
    println(f)

    # Write cell types
    println(f, "CELL_TYPES ", nelems)
    for elem in dom.elems
        println(f, get_vtk_type(elem.shape))
    end
    # Ips cell types
    vtk_vertex = 1
    for ip in ips
        print(f, vtk_vertex, " ")
    end
    println(f)


    if has_data
        # Get node and elem values
        node_vals, node_labels, elem_vals, elem_labels = node_and_elem_vals(dom.nodes, dom.elems)
        nncomps = length(node_labels)
        necomps = length(elem_labels)
    else
        close(f)
        return
    end

    point_data = DTable()
    cell_data  = DTable()

    # Write point data
    println(f, "POINT_DATA ", nnodes)

    # Write vectors
    if :ux in keys(dom.nodes[1].dofdict)
        println(f, "VECTORS ", "|u| float64")
        for node in dom.nodes
            @printf f "%23.10e"   node.dofdict[:ux].U
            @printf f "%23.10e"   node.dofdict[:uy].U
            @printf f "%23.10e\n" dom.ndim==3?node.dofdict[:uz].U:0.0
        end
        for ip in ips
            print(f, "0.0 0.0 0.0\n")
        end
    end
    println(f, )

    # Write nodal scalar data
    for i=1:nncomps
        println(f, "SCALARS ", node_labels[i], " float64 1")
        println(f, "LOOKUP_TABLE default")
        for j in 1:nnodes
            @printf f "%23.10e"  node_vals[j,i]
        end
        for ip in ips
            print(f, "0.0 ")
        end
        println(f, )
    end

    # Get values at ips
    table = DTable()
    for ip in ips
        push!(table, getvals(ip))
    end

    # Write ip scalar data TODO
    for field in table.header
        println(f, "SCALARS ", field, " float64 1")
        println(f, "LOOKUP_TABLE default")
        for elem in dom.elems
            print(f, "0.0 ")
        end
        field_data = table.dict[field]
        for i in 1:nips
            @printf f "%23.10e" field_data[i]
        end
        println(f, )
    end


    # Write element data
    println(f, "CELL_DATA ", ncells)
    for i=1:necomps
        println(f, "SCALARS ", elem_labels[i], " float64 1")
        println(f, "LOOKUP_TABLE default")
        for j in 1:nelems  #naelems
            e_idx = dom.elems[j].id  #aelems
            @printf f "%23.10e"  elem_vals[e_idx,i]
        end
        for ip in ips
            print(f, "0.0 ")
        end
        println(f, )
    end

    # Write elem type
    println(f, "SCALARS cell_type int 1")
    println(f, "LOOKUP_TABLE default")
    for elem in dom.elems
        println(f, elem.shape, " ")
    end
    println(f, )

    close(f)

    if verbose
        println(GREEN, "  file $filename written (Domain)", DEFAULT)
    end

end


function save(nodes::Array{Node,1}, filename::String; dir::Symbol=:nodir, rev::Bool=false, verbose=true)
    # sort nodes
    if dir in (:x, :y, :z)
        nodes = sort(nodes, dir=dir, rev=rev)
    end

    # filling table
    table = DTable()
    dist  = 0.0
    X0    = nodes[1].X
    for node in nodes
        dist += norm(node.X - X0)
        X0    = node.X
        vals  = getvals(node)
        vals[:dist] = dist
        vals[:x]    = node.X[1]
        vals[:y]    = node.X[2]
        vals[:z]    = node.X[3]
        push!(table, vals)
    end

    save(table, filename, verbose)
end

function save(ips::Array{Ip,1}, filename::String; offset::Float64=0.0, dir::Symbol=:nodir, rev::Bool=false, verbose=true)
    # sort ips
    if dir in (:x, :y, :z)
        ips = sort(ips, dir=dir, rev=rev)
    end

    basename, ext = splitext(filename)
    format = (ext == "")? "dat" : ext[2:end]

    if format=="dat"

        # filling table
        table = DTable()
        dist  = offset
        X0    = ips[1].X
        for node in ips
            dist += norm(node.X - X0)
            X0    = node.X
            vals  = getvals(ip.data)
            vals[:dist] = dist
            vals[:x]    = ip.X[1]
            vals[:y]    = ip.X[2]
            vals[:z]    = ip.X[3]
            push!(table, vals)
        end

        save(table, filename, verbose)
        return
    end

    if format=="vtk"
        # Saves the integration points in vtk format

        # Find nips
        nips = length(ips)

        # Open filename
        f = open(filename, "w")

        println(f, "# vtk DataFile Version 3.0")
        println(f, "pyfem output ")
        println(f, "ASCII")
        println(f, "DATASET UNSTRUCTURED_GRID")
        println(f, "")
        println(f, "POINTS ", nips, " float64")

        # Write ip points
        for ip in ips
            @printf f "%23.15e %23.15e %23.15e \n" ip.X[1] ip.X[2] ip.X[3]
        end
        println(f)

        # Write connectivities
        println(f, "CELLS $nips  $(nips*2)")
        for (i,ip) in enumerate(ips)
            println(f, "1 $(i-1)", )
        end
        println(f)

        # Write cell types
        println(f, "CELL_TYPES ", nips)
        vtk_vertex = 1
        for ip in ips
            print(f, vtk_vertex, " ")
        end
        println(f)

        # Get values at ips
        table = DTable()
        for ip in ips
            push!(table, getvals(ip))
        end

        # Write point data
        println(f, "POINT_DATA ", nips)

        # Write ip scalar data
        for field in table.header
            println(f, "SCALARS ", field, " float64 1")
            println(f, "LOOKUP_TABLE default")
            field_data = table.dict[field]
            for i in 1:nips
                @printf f "%23.10e" field_data[i]
            end
            println(f, )
        end

        # Write cell data
        println(f, "CELL_DATA ", nips)

        # Write owner elem type
        println(f, "SCALARS cell_shape_tag int 1")
        println(f, "LOOKUP_TABLE default")
        for ip in ips
            print(f, get_shape_tag(ip.owner.shape), " ")
        end
        println(f, )

        close(f)

        if verbose
            println("  file $filename written (Domain)")
        end
        
    end
end
