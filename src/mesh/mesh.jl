export Mesh
export copy, move
export generate_mesh, save

include("entities.jl")


### Type Mesh
type Mesh
    points ::Array{Point,1}
    cells  ::Array{Cell,1}
    faces  ::Array{Cell,1}
    bpoints::Dict{Uint64,Point}
    ndim   ::Int
    bins   ::Bins

    function Mesh()
        this = new()
        this.points = []
        this.bpoints = Dict{Uint64, Point}()
        this.cells  = []
        this.faces  = []
        this.ndim   = 0
        this.bins   = Bins()
        return this
    end
end


include("block.jl")


function get_surface(cells::Array{Cell,1})
    # Needs cells id to be numbered
    surf_dict = Dict{Uint64, Cell}()

    # Get only unique faces
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

    return [face for (key, face) in surf_dict]
end

function generate_mesh(blocks::Block...; verbose::Bool=true)
    generate_mesh([blocks...], verbose)
end

function generate_mesh(blocks::Array, verbose::Bool=true)
    nblocks = length(blocks)
    if verbose
        println("Mesh generation:")
        println("  analyzing $nblocks block(s)") 
    end
    mesh = Mesh()

    # Split blocks
    for (i,b) in enumerate(blocks)
        b.id = i
        split_block(b, mesh)
        if verbose; print("  spliting block ",i,"...\r") end
    end


    # Get ndim
    ndim = 2
    for point in mesh.points
        if point.z>0.0; ndim = 3; break end
    end
    mesh.ndim = ndim
    
    # Numberig
    for (i,p) in enumerate(mesh.points) p.id = i end
    for (i,c) in enumerate(mesh.cells ) c.id = i; c.ndim=ndim end
    for (i,f) in enumerate(mesh.faces ) f.id = i; f.ndim=ndim end

    # Facets
    if verbose; print("  finding faces...\r") end
    mesh.faces = get_surface(mesh.cells)

    if verbose
        npoints = length(mesh.points)
        ncells  = length(mesh.cells)
        nfaces  = length(mesh.faces)
        println("  $ndim","d found             ")
        @printf "  %5d points obtained\n" npoints
        @printf "  %5d cells obtained\n" ncells
        @printf "  %5d faces obtained\n" nfaces
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


import .Definitions.save
function Definitions.save(mesh::Mesh, filename::String)
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

end


