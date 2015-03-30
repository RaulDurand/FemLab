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

export Mesh
export copy, move
export generate_mesh, save, loadmesh

include("entities.jl")
include("delaunay.jl")

Pkg.installed("JSON")==nothing?  Pkg.add("JSON") : nothing
using JSON


### Type Mesh
type Mesh
    points ::Array{Point,1}
    cells  ::Array{Cell,1}
    faces  ::Array{Cell,1}
    edges  ::Array{Cell,1}
    bpoints::Dict{Uint64,Point}
    ndim   ::Int
    bins   ::Bins

    function Mesh()
        this = new()
        this.points = []
        this.bpoints = Dict{Uint64, Point}()
        this.cells  = []
        this.faces  = []
        this.edges  = []
        this.ndim   = 0
        this.bins   = Bins()
        return this
    end
end


include("block.jl")


function get_surface(cells::Array{Cell,1})
    surf_dict = Dict{Uint64, Cell}()

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in get_faces(cell)
            #conns = [ pt.id for pt in face.points ]
            #hs = hash(conns)
            hs = hash(face)
            if haskey(surf_dict, hs)
                delete!(surf_dict, hs)
            else
                surf_dict[hs] = face
            end
        end
    end

    return [ face for face in values(surf_dict) ]
end

function get_edges(surf_cells::Array{Cell,1})
    edges_dict = Dict{Uint64, Cell}()

    # Get only unique edges
    for cell in surf_cells
        for edge in get_faces(cell)
            edge.ocell = cell.ocell
            edges_dict[hash(edge)] = edge
        end
    end

    return [ edge for edge in values(edges_dict) ] 
end

function get_surface_alt(cells::Array{Cell,1})
    # Actually slower....
    # Get all points
    pointsd = Dict{Uint64, Point}()
    for cell in cells
        for point in cell.points
            pointsd[hash(point)] = point
        end
    end
    points = values(pointsd)

    # Get incidence matrix (shares) (fast)
    np = length(points)
    N = [ Cell[] for i=1:np]
    for cell in cells
        for pt in cell.points
            push!(N[pt.id], cell)
        end
    end

    @show 1

    # Get matrix of cells faces
    F = [ get_faces(cell) for cell in cells]
    nc = length(cells)
    #CF = Array(Array{Array{Int64,1},1}, nc)
    CF = Array(Array{Uint64,1}, nc)
    for cell in cells # fast
        #CF[cell.id] = [ sort([pt.id for pt in face.points]) for face in F[cell.id]]
        CF[cell.id] = [ hash(face) for face in F[cell.id]]
    end
    @show 2

    # Get cells boundary flag matrix
    CB = [ trues(length(CF[cell.id])) for cell in cells]
    for cell in cells
        for (i,fcon) in enumerate(CF[cell.id])
            for pid in fcon
                for cl in N[pid]
                    if cl.id == cell.id; continue end
                    if fcon in CF[cl.id]
                        CB[cell.id][i] = false
                    end
                end
            end
        end
    end
    @show 3

    # Get list of boundary faces (almost fast)
    facets = Cell[]
    for cell in cells
        for (i,face) in enumerate(F[cell.id])
            if CB[cell.id][i]
                push!(facets, face)
            end
        end
    end
    @show 4
    return facets
end

function generate_mesh(blocks::Block...; verbose::Bool=true, genfacets=true, genedges=false, initial_mesh=nothing)
    generate_mesh([blocks...], verbose, genfacets, genedges, initial_mesh)
end

function generate_mesh(initial_mesh::Mesh, blocks::Block...; verbose::Bool=true, genfacets=true, genedges=false)
    generate_mesh([blocks...], verbose, genfacets, genedges, initial_mesh)
end

function generate_mesh(blocks::Array, verbose::Bool=true, genfacets=true, genedges=false, initial_mesh=nothing)
    nblocks = length(blocks)
    if verbose
        println(BOLD, CYAN, "Mesh generation:", DEFAULT)
        println("  analyzing $nblocks block(s)") 
    end

    mesh = initial_mesh==nothing? Mesh(): initial_mesh

    # Split blocks
    for (i,b) in enumerate(blocks)
        b.id = i
        split_block(b, mesh)
        if verbose; print("  spliting block ",i,"...\r") end
    end

    # Get ndim
    ndim = 2
    for point in mesh.points
        if point.z != 0.0; ndim = 3; break end
    end
    mesh.ndim = ndim
    
    # Numberig
    for (i,p) in enumerate(mesh.points) p.id = i end
    for (i,c) in enumerate(mesh.cells ) c.id = i; c.ndim=ndim end

    # Facets
    if genfacets
        if verbose; print("  finding faces...\r") end
        mesh.faces = get_surface(mesh.cells)
    end

    # Edges
    if genedges && ndim==3
        if verbose; print("  finding edges...\r") end
        mesh.edges = get_edges(mesh.faces)
    end

    if verbose
        npoints = length(mesh.points)
        ncells  = length(mesh.cells)
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        println("  $ndim","d found             ")
        @printf "  %5d points obtained\n" npoints
        @printf "  %5d cells obtained\n" ncells
        if genfacets
            @printf "  %5d faces obtained\n" nfaces
        end
        if genedges
            @printf "  %5d surface edges obtained\n" nedges
        end
        icount = 0
        for block in blocks
            if typeof(block) == BlockInset; icount += block.icount end
        end
        if icount>0
            @printf "  %5d intersections found\n" icount
        end
    end

    return mesh
end


#import .Definitions.save
function save(mesh::Mesh, filename::String, verbose=true)
    # Saves the mesh information in vtk format

    npoints = length(mesh.points)
    ncells  = length(mesh.cells)

    ndata = 0
    for cell in mesh.cells
        ndata += 1 + length(cell.points)
    end

    has_crossed = any([cell.crossed for cell in mesh.cells])

    f = open(filename, "w")

    println(f, "# vtk DataFile Version 3.0")
    println(f, "pyfem output ")
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", npoints, " float64")

    # Write points
    for point in mesh.points
        @printf f "%23.15e %23.15e %23.15e \n" point.x point.y point.z
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ",ncells, " ", ndata)
    for cell in mesh.cells
        print(f, length(cell.points), " ")
        for point in cell.points
            print(f, point.id-1, " ")
        end
        println(f)
    end
    println(f)

    # Write cell types
    println(f, "CELL_TYPES ", ncells)
    for cell in mesh.cells
        println(f, get_vtk_type(cell.shape))
    end
    println(f)

    println(f, "CELL_DATA ", ncells)

    # Write cell type as cell data
    println(f, "SCALARS cell_type int 1")
    println(f, "LOOKUP_TABLE default")
    for cell in mesh.cells
        println(f, cell.shape)
    end
    println(f, )

    # Write flag for crossed cells
    if has_crossed
        println(f, "SCALARS crossed int 1")
        println(f, "LOOKUP_TABLE default")
        for cell in mesh.cells
            println(f, int(cell.crossed))
        end
        println(f, )
    end

    if verbose
        println(GREEN, "  file $filename written (Mesh)", DEFAULT)
    end

end

function get_shape_type(geo, npoints=0)
    types = [ LIN2, LIN3, -1, TRI3, TRI6, -1, QUAD4, QUAD8, QUAD9, TET4, TET10, HEX8, HEX20, -2, LIN4, TRI10, QUAD12, QUAD16]
    shapetype = types[geo+1]
    if shapetype==-2 #joint
        shapetype = npoints==2? LINK2: LINK3
    end

    return shapetype
end


function loadmesh(filename; format="json")
    mesh = Mesh()

    data = JSON.parsefile(filename)
    verts = data["verts"]
    cells = data["cells"]

    # Loading points
    for vert in verts
        C = float(vert["c"])
        C = [C, 0.0][1:3]
        #push!(C, 0.0)
        P = Point(C)
        P.id  = vert["id"]+1
        P.tag = string(vert["tag"])
        push!(mesh.points, P)

    end

    # Get ndim
    ndim = 2
    for point in mesh.points
        if point.z != 0.0; ndim = 3; break end
    end
    mesh.ndim = ndim

    # Loading cells
    for cell in cells
        geo   = cell["geo"]
        conn  = cell["verts"]
        conn  = [ i+1 for i in conn]
        npoints = length(conn)

        # check for embedded joint
        if haskey(cell, "jlinid") # is joint element
            lincell = cells[cell["jlinid"]]
            npoints = length(lincell["verts"])
        end

        shapetype = get_shape_type(geo, npoints)
        points    = [ mesh.points[i] for i in conn ]
        C = Cell(shapetype, points, string(cell["tag"]))
        C.id  = cell["id"]+1
        push!(mesh.cells, C)
    end

    # Generating faces
    mesh.faces = get_surface(mesh.cells)
    surf_dict = [ hash(F) => F for F in mesh.faces]

    all_faces = Face[]
    for (i,C) in enumerate(mesh.cells)
        if is_solid(C.shape)
            faces = get_faces(C)
            ftags = get(cells[i],"ftags", [])
            #@show ftags
            if length(ftags) > 0
                for (j,F) in enumerate(faces)
                    tag = string(ftags[j])
                    if ftags != "0"
                        F.tag = tag
                        surf_dict[hash(F)] = F
                    end
                end
            end
        end
    end

    mesh.faces = [ face for face in values(surf_dict) ]

    return mesh
end
