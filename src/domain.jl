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


"""
`Domain(mesh, [filekey="out"])`

Creates an `Domain` object based on a Mesh object `mesh` and represents the geometric domain to be analysed by the finite element analysis.

**Fields**

`nodes`: An array of nodes

`elems`: An array of finite elements

`faces`: An array of `Face` objects containing the boundary faces

`edges`: An array of `Edge` objects containing all the boundary faces

`node_bcs`: An array of `NodeBC` objects containing all nodal boundary conditions

`face_bcs`: An array of `FaceBC` objects conditions all face boundary conditions

`edge_bcs`: An array of `EdgeBC` objects conditions all edge boundary conditions

`filekey` : An string object that is used as part of the filename of resulting analyses files

"""
mutable struct Domain
    ndim ::Int
    nodes::Array{Node,1}
    elems::Array{Element,1}
    faces::Array{Face,1}
    edges::Array{Edge,1}
    bcs::Array{BC,1}
    filekey::AbstractString

    monitors::Array{Monitor,1}
    nincs::Integer
    nouts::Integer
    ndofs::Integer

    function Domain(;filekey::AbstractString="out")
        this = new()

        this.bcs      = []
        this.filekey  = filekey

        this.monitors = []
        this.nincs    = 0
        this.ndofs    = 0
        this.nouts    = 0
        return this
    end
end


function Domain(mesh::Mesh; filekey::AbstractString="out", stress_state=:general, thickness=1.0)
    dom  = Domain(filekey=filekey)
    ndim = dom.ndim = mesh.ndim
    global gl_stress_state = stress_state
    global gl_thickness    = thickness

    # Setting nodes
    dom.nodes = [ Node(point, id=i) for (i,point) in enumerate(mesh.points)]

    # Setting elements
    dom.elems = Array{Element}(0)
    for (i,cell) in enumerate(mesh.cells)
        conn = [ p.id for p in cell.points ]
        elem = Element(cell.shape, dom.nodes[conn], ndim, cell.tag)
        elem.id = i
        push!(dom.elems, elem)
    end

    # Setting linked cells
    for (i,cell) in enumerate(mesh.cells)
        # setting linked elements
        for lcell in cell.linked_cells
            id = lcell.id
            push!(dom.elems[i].linked_elems, dom.elems[id])
        end
    end

    # Setting faces
    dom.faces = Array{Face}(0)
    for cell in mesh.faces
        conn = [ p.id for p in cell.points ]
        face = Face(cell.shape, dom.nodes[conn], ndim)
        face.oelem = dom.elems[cell.ocell.id]
        push!(dom.faces, face)
    end

    # Setting edges
    dom.edges = Array{Edge}(0)
    for cell in mesh.edges
        conn = [ p.id for p in cell.points ]
        edge = Edge(cell.shape, dom.nodes[conn], ndim)
        edge.oelem = dom.elems[cell.ocell.id]
        push!(dom.edges, edge)
    end

    return dom
end


function dom_load_json(filename::AbstractString)
    mesh = Mesh(filename)
    dom  = Domain(mesh)

    file = open(filename)
    data = JSON.parse(readall(file))
    close(file)

    # Read materials list
    mats  = Material[]
    for dmat in data["materials"] # array of materials data
        id    = dmat["id"]
        model = eval(parse(dmat["model"]))
        delete!(dmat, "id")
        delete!(dmat, "model")

        args = Dict( Symbol(k) => v for (k,v) in dmat )
        mat = model(;args...)
        push!(mats, mat)
    end

    # Read and set material for each element
    cells = data["cells"]
    for i in 1:endof(cells)
        mat_id = cells[i]["mat"]
        set_mat( dom.elems[i], mats[mat_id] )
    end

    # Node boundary conditions
    node_bc = data["node_bc"]
    for bc in node_bc
        id = get(bc,"id",0)

        if id>0
            delete!(bc, "id")
            args = Dict( Symbol(k) => v for (k,v) in  bc )
            apply_bc(dom.nodes[id]; args...)
        else
            cond = parse(bc["cond"])
            delete!(bc, "cond")
            args = Dict( Symbol(k) => v for (k,v) in  bc )
            apply_bc( dom.nodes[cond]; args...)
        end
    end

    # Face boundary conditions
    face_bc = data["face_bc"]
    for bc in face_bc
        if haskey(bc, "id")
            error("Unsupported key id for face boundary condition while reading input file.")
        end

        cond = parse(bc["cond"])
        delete!(bc, "cond")
        args = Dict( Symbol(k) => v for (k,v) in  bc )
        apply_bc( dom.faces[cond]; args...)
    end

    return dom
end


function Domain(filename::AbstractString; filekey::AbstractString="out")
    ext  = splitext(filename)[2]

    if ext in (".json", ".msh")
        dom = dom_load_json(filename)
        dom.filekey = filekey
        return dom
    else
        mesh = Mesh(filename)
        return Domain(mesh, filekey=filekey)
    end

end


function set_bc(dom::Domain, bcs::BC...)
    dom.bcs = []
    for bc in bcs
        setup_bc!(bc, dom)
        push!(dom.bcs, bc)
    end
end

function set_monitors(dom::Domain, monitors::Monitor...)
    dom.monitors = []
    for mon in monitors
        setup_monitor!(mon, dom)
        push!(dom.monitors, mon)
    end
end
set_trackers = set_monitors

function update_monitors!(dom::Domain)
    for monitor in dom.monitors
        update_monitor!(monitor)
    end
end


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
        #patch = node.shares[ [issolid(elem) for elem in node.shares]  ]
        #patch = node.shares[ map(is_solid, node.shares) ]
    #end
end


function state_backup(elems::Array{Element,1})
    for elem in elems
        for ip in elem.ips
            ip.data_bk = deepcopy(ip.data)
        end
    end
end

function state_restore(elems::Array{Element,1})
    for elem in elems
        for ip in elem.ips
            ip.data = deepcopy(ip.data_bk)
        end
    end
end


function all_node_and_elem_vals(nodes::Array{Node,1}, elems::Array{Element,1})

    # Calculate max number of node and elem labes
    matset  = Set{Material}()
    nlabels = Set{Symbol}()
    elabels = Set{Symbol}()
    for elem in elems
        if !(elem.mat in matset)
            # getting values from current element
            node_vals, elem_vals = elem_and_node_vals(elem.mat, elem)
            union!(nlabels, keys(node_vals))
            union!(elabels, keys(elem_vals))
            push!(matset, elem.mat)
        end
    end

    nlabels_idx = Dict( key=>i for (i,key) in enumerate(nlabels) )
    elabels_idx = Dict( key=>i for (i,key) in enumerate(elabels) )
    nncomps     = length(nlabels)  # number of nodal labels
    necomps     = length(elabels)  # number of elem labels

    nnodes     = length(nodes)
    nelems     = length(elems)
    Node_vals  = zeros(Float32, nnodes, nncomps) # float32 is important for paraview
    Node_reps  = zeros(Int64  , nnodes, nncomps)
    Elem_vals  = zeros(Float32, nelems, necomps) # float32 is important for paraview

    # In elements
    for elem in elems
        # getting values from current element
        node_vals, elem_vals = elem_and_node_vals(elem.mat, elem)

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
            node_vals, elem_vals = elem_and_node_vals(elem.mat, elem)
            union!(nlabels, keys(node_vals))
            union!(elabels, keys(elem_vals))
            push!(matset, elem.mat)
        end
    end

    nlabels_idx = Dict( key=>i for (i,key) in enumerate(nlabels) )
    elabels_idx = Dict( key=>i for (i,key) in enumerate(elabels) )
    nncomps     = length(nlabels)
    necomps     = length(elabels)

    nnodes     = length(nodes)
    nelems     = length(elems)
    Node_vals  = zeros(Float32, nnodes, nncomps) # float32 is important for paraview
    Node_reps  = zeros(Int64  , nnodes, nncomps)

    # In elements
    for elem in elems
        # getting values from current element
        node_vals, elem_vals = elem_and_node_vals(elem.mat, elem)

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


function save(dom::Domain, filename::AbstractString; verbose=true, save_ips=false)
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
    println(f, "FemLab output ")
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
        println(f, elem.shape.vtk_type)
    end
    println(f)

    # check if all elements have material defined
    has_data = all([ isdefined(elem, :mat) for elem in dom.elems])

    if has_data
        # Get node and elem values
        node_vals, node_labels, elem_vals, elem_labels = all_node_and_elem_vals(dom.nodes, dom.elems)
        nncomps = length(node_labels)
        necomps = length(elem_labels)
    else
        close(f)
        verbose && print_with_color(:green, "  file $filename written (Domain)\n")
        return
    end

    # Write point data
    println(f, "POINT_DATA ", nnodes)

    # Write vectors
    if :uy in node_labels
        ux_idx = findfirst(node_labels, :ux)
        uy_idx = findfirst(node_labels, :uy)
        uz_idx = findfirst(node_labels, :uz)
        println(f, "VECTORS ", "U float64")
        for (i,node) in enumerate(dom.nodes)
            @printf f "%23.10e"   node_vals[i,ux_idx]
            @printf f "%23.10e"   node_vals[i,uy_idx]
            @printf f "%23.10e\n" uz_idx==0?0.0:node_vals[i,uz_idx]
        end
    end
    println(f, )

    # Write node numbers
    println(f, "SCALARS ", "Node-ID", " int 1")
    println(f, "LOOKUP_TABLE default")
    for node in dom.nodes
        @printf f "%5d " node.id
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

    # Write elem numbers
    println(f, "SCALARS ", "Elem-ID", " int 1")
    println(f, "LOOKUP_TABLE default")
    for elem in dom.elems  #naelems
        @printf f "%5d " elem.id
    end
    println(f, )

    # Write element scalar data
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
        println(f, elem.shape.vtk_type, " ")
    end
    println(f, )

    close(f)

    verbose && print_with_color(:green, "  file $filename written (Domain)\n")

    # save ip information as vtk
    if has_data && save_ips
        basename, ext = splitext(filename)
        save_dom_ips(dom, basename*"-ip"*ext, verbose)
    end
end


# Saves ips information from domain as a vtk file
function save_dom_ips(dom::Domain, filename::AbstractString, verbose=true)

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
        push!(table, ip_vals(ip))
    end

    # Write point data
    println(f, "POINT_DATA ", nips)

    # Write ip scalar data TODO
    for field in keys(table)
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
        print_with_color(:green, "  file $filename written (Domain)\n")
    end
end


function save2(dom::Domain, filename::AbstractString; verbose=true)
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
        println(f, elem.shape.vtk_type)
    end
    # Ips cell types
    vtk_vertex = 1
    for ip in ips
        print(f, vtk_vertex, " ")
    end
    println(f)

    if has_data
        # Get node and elem values
        node_vals, node_labels, elem_vals, elem_labels = all_node_and_elem_vals(dom.nodes, dom.elems)
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
        push!(table, ip_vals(ip))
    end

    # Write ip scalar data TODO
    for field in keys(table)
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
        print_with_color(:green, "  file $filename written (Domain)\n")
    end
end


function save(nodes::Array{Node,1}, filename::AbstractString; dir::Symbol=:nodir, rev::Bool=false, verbose::Bool=true)
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
        vals  = node_vals(node)
        vals[:dist] = dist
        vals[:x]    = node.X[1]
        vals[:y]    = node.X[2]
        vals[:z]    = node.X[3]
        push!(table, vals)
    end

    save(table, filename, verbose)
end


function save(ips::Array{Ip,1}, filename::AbstractString; offset::Float64=0.0, dir::Symbol=:nodir, rev::Bool=false, verbose=true)
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
            vals  = ip_state_vals(ip.data)
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
            push!(table, ip_vals(ip))
        end

        # Write point data
        println(f, "POINT_DATA ", nips)

        # Write ip scalar data
        for field in keys(table)
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
            print_with_color(:green, "  file $filename written (Domain Ips)\n")
        end
    end
end
