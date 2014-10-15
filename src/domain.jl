
#using  Entities
#include("entities.jl")
#using  .MeshGen
export Domain
export track

type Face
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union(Element,Nothing)
    function Face(shape, nodes, ndim)
        this = new(shape, nodes, ndim)
        this.oelem = nothing
        return this
    end
end

#function getindex(faces::Array{Face,1}, cond::Expr)
    #result = Array(Face,0)
    #for face in faces
        #mcond  = copy(cond)
        #coords = getcoords(face.nodes)
        #x = all(coords[:,1].==coords[1,1]) ? coords[1,1] : NaN
        #y = all(coords[:,2].==coords[1,2]) ? coords[1,2] : NaN
        #z = all(coords[:,3].==coords[1,3]) ? coords[1,3] : NaN
        #substitute!(mcond, (:x, :y, :z), (x, y, z))
        #if eval(mcond); push!(result, face) end
    #end
    #return result
#end

function getindex(faces::Array{Face,1}, cond::Expr)
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = cond
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Node getindex: Invalid condition ", cond)
    end

    result = Array(Face,0)
    for face in faces
        coords = getcoords(face.nodes)
        x = all(coords[:,1].==coords[1,1]) ? coords[1,1] : NaN
        y = all(coords[:,2].==coords[1,2]) ? coords[1,2] : NaN
        z = all(coords[:,3].==coords[1,3]) ? coords[1,3] : NaN
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
    if length(faces)==0; println("Warning, applying boundary conditions to empty array of faces") end
    for face in faces
        set_bc(face; args...)
    end
end


typealias Edge Face

#type TrackStruct{T}
    #obj::T
    #data::DTable
    #function TrackStruct(obj::T, data::Dtable)
        #this = new(obj, data)
        #return this
    #end
#end

type Domain
    ndim ::Int
    nodes::Array{Node,1}
    elems::Array{Element,1}
    faces::Array{Face,1}
    edges::Array{Edge,1}
    trk_nodes::Array{(Node, DTable), 1}
    trk_ips  ::Array{(Ip, DTable), 1}
    function Domain(mesh::Mesh)
        this = new()
        load_mesh(this, mesh)
        this.trk_nodes = []
        this.trk_ips   = []
        return this
    end
end

function get_faces(elem::Element)
    faces  = Array(Face,0)
    if !haskey(FACETS_IDXS, elem.shape)
        return faces
    end

    f_idxs = FACETS_IDXS[elem.shape]

    for f_idx in f_idxs
        points = [ elem.points[i] for i in f_idx]
        face   = Cell(FACETS_SHAPE[elem.shape], points)
        push!(faces, face)
    end

    return faces
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
    #@show length(dom.elems)

    # Setting faces
    dom.faces = Array(Face,0)
    for cell in mesh.faces
        conn = [ p.id for p in cell.points ]
        face = Face(cell.shape, dom.nodes[conn], ndim)
        #@show cell.id
        face.oelem = dom.elems[cell.ocell.id]
        push!(dom.faces, face)
    end

    # Setting edges
    dom.edges = [] #TODO

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

#function track(dom::Domain, node::Node)
    #headr = vcat([ [dof.sU, dof.sF] for dof in node.dofs]...)
    #headr = [:id, :time, headr] 
    #table = DTable(headr)
    #table.extra[:entity] = node
    #push!(dom.track_tables, table)
    #return table
#end

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


function tracking(dom::Domain)
    # tracking ips
    for (ip, table) in dom.trk_ips
        vals = getvals(ip.data)
        push!(table, vals)
    end

    # tracking nodes
    for (node, table) in dom.trk_nodes
        vals = getvals(node)
        push!(table, vals)
    end

end


import .Definitions.save
function save(dom::Domain, filename::String)
    # Saves the dom information in vtk format

    nnodes = length(dom.nodes)
    nelems  = length(dom.elems)

    # Number of total connectivities
    nconns = 0
    for elem in dom.elems
        nconns += 1 + length(elem.nodes)
    end

    # Get node and elem values
    node_vals, node_labels, elem_vals, elem_labels = node_and_elem_vals(dom.nodes, dom.elems)
    nncomps = length(node_labels)
    necomps = length(elem_labels)

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

    # Write cell tag
    #println("SCALARS tag int 1")
    #println("LOOKUP_TABLE default")
    #for j=1:naelems
        #try
            #tag = int(dom.aelems[j].tag)
        #except ValueError
            #tag = 0
        #println("%d"%(tag))
    #println()

    #if not dom.track_per_inc
        #dom.write_history(node_vals, node_labels)

end


function save(nodes::Array{Node,1}, filename::String; dir::Symbol=:nodir, rev::Bool=false)
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

    save(table, filename)
end

function save(ips::Array{Ip,1}, filename::String; offset::Float64=0.0, dir::Symbol=:nodir, rev::Bool=false)
    # sort ips
    if dir in (:x, :y, :z)
        ips = sort(ips, dir=dir, rev=rev)
    end

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

    save(table, filename)
end
