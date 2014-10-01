export Mesh, Block2D, Block3D, BlockTruss, BlockInset
export generate_mesh, save

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


### Type Block

abstract Block

type Block2D <: Block
    coords::Array{Float64,2}
    nx::Int64
    ny::Int64
    shape::ShapeType
    tag::String
    id::Int64
    function Block2D(coords; nx=1, ny=1, shape=QUAD4, tag="", id=-1)
        C    = size(coords,1)==2? box_coords(vec(coords[1,:]), vec(coords[2,:])): coords
        this = new(C, nx, ny, shape, tag, id)
        return this
    end
end

type Block3D <: Block
    coords::Array{Float64,2}
    nx::Int64
    ny::Int64
    nz::Int64
    shape::ShapeType
    tag::String
    id::Int64
    function Block3D(coords; nx=1, ny=1, nz=1, shape=HEX8, tag="", id=-1)
        C    = size(coords,1)==2? box_coords(vec(coords[1,:]), vec(coords[2,:])): coords
        this = new(C, nx, ny, nz, shape, tag, id)
        return this
    end
end

#type BlockInset <: Block
    #coords::Array{Float64,2}
    #curvetype::Int64
    #closed::Bool
    #shape::ShapeType
    #tag::String
    #id::Int64
    #function BlockInset(coords; curvetype=3, closed=false, shape=LIN3, tag="", id=-1)
        #this = new(coords, curvetype, closed, shape, tag, id)
        #return this
    #end
#end


function box_coords{T1<:Number, T2<:Number}(C1::Array{T1,1}, C2::Array{T2,1})
    C = Array(Float64, 8, 3)
    x1 = C1[1]
    y1 = C1[2]
    lx = C2[1] - C1[1]
    ly = C2[2] - C1[2]

    if length(C1)==2
        return [
            x1     y1     0.0
            x1+lx  y1     0.0
            x1+lx  y1+ly  0.0
            x1     y1+ly  0.0 ]
    else
        z1 = C1[3]
        lz = C2[3] - C1[3]
        return [
            x1      y1      z1 
            x1+lx   y1      z1 
            x1+lx   y1+ly   z1 
            x1      y1+ly   z1 
            x1      y1      z1+lz 
            x1+lx   y1      z1+lz 
            x1+lx   y1+ly   z1+lz 
            x1      y1+ly   z1+lz ]
    end
end


function split_block(bl::Block2D, msh::Mesh)
    nx, ny = bl.nx, bl.ny
    shape = bl.shape

    if shape==QUAD4
        p_arr = Array(Point, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = shape_func(shape, [r, s])
                C = round(N'*bl.coords, 8)
                C = reshape(C, 3)
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    p = get_point(msh.bpoints, C)
                    if p==nothing
                        p = Point(C); push!(msh.points, p)
                        msh.bpoints[hash(p)] = p
                    end
                else
                    p = Point(C); push!(msh.points, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:ny
            for i = 1:nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+1, j  ]
                p3 = p_arr[i+1, j+1]
                p4 = p_arr[i  , j+1]

                cell = Cell(shape, [p1, p2, p3, p4], bl.tag)
                push!(msh.cells, cell)
                # Faces
                #faces_con = Array(Point, 0)
                #if i==1;  push!(faces_con, [p3, p0]); end
                #if i==nx; push!(faces_con, [p1, p2]); end
                #if j==1;  push!(faces_con, [p0, p1]); end
                #if j==ny; push!(faces_con, [p2, p3]); end
                #for fc in faces_con
                    #face = Cell(LIN2, fc)
                    #push!(msh.faces, face)
                #end
            end
        end
    elseif shape == QUAD8 || shape == QUAD9
        p_arr = Array(Point, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                if shape==QUAD8 && iseven(i) && iseven(j) continue end

                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = shape_func(shape, [r, s])
                C = round(N'*bl.coords, 8)
                C = reshape(C, 3)
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    p = get_point(msh.bpoints, C)
                    if p==nothing
                        p = Point(C); push!(msh.points, p)
                        msh.bpoints[hash(p)] = p
                    end
                else
                    p = Point(C); push!(msh.points, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:2:2*ny+1
            for i = 1:2:2*nx+1
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p4 = p_arr[i  , j+2]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p7 = p_arr[i+1, j+2]
                p8 = p_arr[i  , j+1]

                if shape==QUAD8
                    cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8], bl.tag)
                else
                    p9   = p_arr[i+1, j+1]
                    cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8, p9], bl.tag)
                end
                push!(msh.cells, cell)
            end
        end

    end
end

function split_block(bl::Block3D, msh::Mesh)
    nx, ny, nz = bl.nx, bl.ny, bl.nz
    shape = bl.shape

    if shape==HEX8
        p_arr = Array(Point, nx+1, ny+1, nz+1)
        for k = 1:nz+1
            for j = 1:ny+1
                for i = 1:nx+1
                    r = (2.0/nx)*(i-1) - 1.0
                    s = (2.0/ny)*(j-1) - 1.0
                    t = (2.0/nz)*(k-1) - 1.0
                    N = shape_func(shape, [r, s, t])
                    C = round(N'*bl.coords, 8)
                    C = reshape(C, 3)
                    p::Any = nothing
                    if i in (1, nx+1) || j in (1, ny+1) || k in (1, nz+1)
                        p = get_point(msh.bpoints, C)
                        if p==nothing
                            p = Point(C); push!(msh.points, p)
                            msh.bpoints[hash(p)] = p
                        end
                    else
                        p = Point(C); push!(msh.points, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    conn = [
                        p_arr[i  , j  , k  ],
                        p_arr[i+1, j  , k  ],
                        p_arr[i+1, j+1, k  ],
                        p_arr[i  , j+1, k  ],
                        p_arr[i  , j  , k+1],
                        p_arr[i+1, j  , k+1],
                        p_arr[i+1, j+1, k+1],
                        p_arr[i  , j+1, k+1]]

                    cell = Cell(shape, conn, bl.tag)
                    push!(msh.cells, cell)
                end
            end
        end
    elseif shape == HEX20
        p_arr = Array(Point, 2*nx+1, 2*ny+1, 2*nz+1)
        for k = 1:2*nz+1
            for j = 1:2*ny+1
                for i = 1:2*nx+1
                    if iseven(i) && iseven(j) continue end
                    if iseven(j) && iseven(k) continue end
                    if iseven(k) && iseven(i) continue end

                    r = (2.0/nx)*(i-1) - 1.0
                    s = (2.0/ny)*(j-1) - 1.0
                    t = (2.0/nz)*(k-1) - 1.0
                    N = shape_func(shape, [r, s])
                    C = round(N'*bl.coords, 8)
                    C = reshape(C, 3)
                    p::Any = nothing
                    if i in (1, 2*nx+1) || j in (1, 2*ny+1) || k in (1, 2*nz+1)
                        p = get_point(msh.bpoints, C)
                        if p==nothing
                            p = Point(C); push!(msh.points, p)
                            msh.bpoints[hash(p)] = p
                        end
                    else
                        p = Point(C); push!(msh.points, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:2:2*nz+1
            for j = 1:2:2*ny+1
                for i = 1:2:2*nx+1
                    conn = [
                        p_arr[i  , j  , k  ],
                        p_arr[i+2, j  , k  ],
                        p_arr[i+2, j+2, k  ],
                        p_arr[i  , j+2, k  ],
                        p_arr[i  , j  , k+2],
                        p_arr[i+2, j  , k+2],
                        p_arr[i+2, j+2, k+2],
                        p_arr[i  , j+2, k+2],
                                            
                        p_arr[i+1, j  , k  ],
                        p_arr[i+2, j+1, k  ],
                        p_arr[i+1, j+2, k  ],
                        p_arr[i  , j+1, k  ],
                        p_arr[i+1, j  , k+2],
                        p_arr[i+2, j+1, k+2],
                        p_arr[i+1, j+2, k+2],
                        p_arr[i  , j+1, k+2],
                                           
                        p_arr[i  , j  , k+1],
                        p_arr[i+2, j  , k+1],
                        p_arr[i+2, j+2, k+1],
                        p_arr[i  , j+2, k+1]]

                    cell = Cell(shape, conn, bl.tag)
                    push!(msh.cells, cell)
                end
            end
        end

    end
end


function get_surface(cells::Array{Cell,1})
    surf_dict = Dict{Uint64, Cell}()

    # Get only unique faces
    for cell in cells
        for face in get_faces(cell)
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
    
    # Facets
    if verbose; print("  finding faces...\r") end
    mesh.faces = get_surface(mesh.cells)

    # Numberig
    for (i,p) in enumerate(mesh.points) p.id = i end
    for (i,c) in enumerate(mesh.cells ) c.id = i; c.ndim=ndim end
    for (i,f) in enumerate(mesh.faces ) f.id = i; f.ndim=ndim end

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



type BlockTruss <: Block
    coords::Array{Float64,2}
    conns ::Array{Int64,2}
    shape ::ShapeType
    tag::String
    id::Int64

    function BlockTruss(coords::Array{Float64,2}, conns::Array{Int64,2}; shape=LIN2, tag="", id=-1)
        ncols = size(coords,2)
        if !(ncols in (2,3)); error("Invalid coordinates matrix for BlockTruss") end
        if ncols==2
            C = [coords  zeros(size(coords,1)) ]
        else
            C = coords
        end
        this = new(C, conns, shape, tag, id)
        return this
    end
end

function split_block(bl::BlockTruss, msh::Mesh)
    n = size(bl.coords, 1) # number of points
    m = size(bl.conns , 1) # number of truss cells
    p_arr = Array(Point, n)
    for i=1:n
        C = reshape(bl.coords[i,:], 3)
        p = get_point(msh.bpoints, C)
        if p==nothing; 
            p = Point(C) 
            msh.bpoints[hash(p)] = p
            push!(msh.points, p)
        end
        p_arr[i] = p
    end
    for i=1:m
        p1 = p_arr[bl.conns[i, 1]]
        p2 = p_arr[bl.conns[i, 2]]
        cell = Cell(bl.shape, [p1, p2], bl.tag)
        push!(msh.cells, cell)
    end
end

type BlockInset <: Block
    coords   ::Array{Float64,2}
    curvetype::Union(Int,String) # 0:polyline, 1:closed polyline, 2: lagrangian, 3:cubic Bezier with inner points
    closed   ::Bool
    shape    ::ShapeType
    tag      ::String
    id       ::Int64
    icount   ::Int64
    _endpoint  ::Union(Point, Nothing)
    _startpoint::Union(Point, Nothing)
    function BlockInset(coords; curvetype=0, closed=false, shape=LIN3, tag="", id=-1)
        if typeof(curvetype)<:Integer
            if !(0<=curvetype<=3); error("Wrong curve type") end
            ctype = curvetype
        else
            cases = ["polyline"=>0, "closed polyline"=>1, "lagrangian"=>2, "Bezier"=>3, "bezier"=>3]
            ctype = get(cases, curvetype, -1)
            if ctype==-1; error("Wrong curve type") end
        end
        this = new(coords, ctype, closed, shape, tag, id)
        this.icount = 0
        this._endpoint   = nothing
        this._startpoint = nothing
        return this
    end
end

#import Definitions.@out

function cubicBezier(s::Float64, PorQ::Array{Float64,2}, isQ=false)
    # check
    if size(PorQ,1) < 4
        println("PorQ = ", PorQ)
        error("List of points must have at least 4 points for cubic Bezier.")
    end

    # input data
    ndim = size(PorQ,2)    # space dimension
    np   = size(PorQ,1)    # number of points
    ns   = np - 1          # number of spacings
    nb   = ifloor(ns/3)    # number of bezier curves
    Ds   = 1.0 / nb        # spacing between curves

    # find index of Bezier and local coordinate t
    ib = ifloor(s/Ds) + 1    # index of Bezier
    if ib > nb; ib = nb end # fix index if s ~= 1+eps
    s0 = (ib-1) * Ds        # s @ left point
    t  = (s - s0) / Ds      # local t
    if t > 1.0; t = 1.0 end # clean rubbish. e.g. 1.000000000000002

    # collect control points
    Q = zeros(4, ndim)    # control points
    k = 1 + (ib-1) * 3    # position of first point of bezier
    if isQ
        for i in 1:4
            Q[i,:] = PorQ[k+i-1]
        end
    else
        PQ1 = PorQ[k  ,:]
        PQ2 = PorQ[k+1,:]
        PQ3 = PorQ[k+2,:]
        PQ4 = PorQ[k+3,:]
        Q[1,:] =         PQ1
        Q[2,:] = (-5.0 * PQ1 + 18.0 * PQ2 -  9.0 * PQ3 + 2.0 * PQ4) / 6.0
        Q[3,:] = ( 2.0 * PQ1 -  9.0 * PQ2 + 18.0 * PQ3 - 5.0 * PQ4) / 6.0
        Q[4,:] =                                               PQ4
    end

    # compute Bezier
    Q1, Q2, Q3, Q4 = Q[1,:], Q[2,:], Q[3,:], Q[4,:]
    a =       Q4 - 3.0 * Q3 + 3.0 * Q2 - Q1
    b = 3.0 * Q3 - 6.0 * Q2 + 3.0 * Q1
    c = 3.0 * Q2 - 3.0 * Q1
    d =       Q1
    return vec(a*t*t*t + b*t*t + c*t + d)
end

function interLagrange(s, coords)
    #  Interpolates coordinates for s between 0 and 1
    #
    #  0                   +1
    #  0---1---2---3---..---n  -->s
    n = size(coords, 1)
    S = zeros(n)

    for i in 1:n
        f  = 1.0
        si = 1.0*(i-1)/(n-1.0)
        for j in 1:n
            sj = 1.0*(j-1)/(n-1.0)
            if i != j
                f *= (s - sj)/(si - sj)
            end
        end
        S[i] = f
    end
    return coords'*S
end


function split_block(bl::BlockInset, msh::Mesh)
    n, ndim = size(bl.coords)

    if n<2; error("At list two points are required in BlockInset") end
    # 0:polyline, 1:closed polyline, 2: lagrangian, 3:cubic Bezier, 4:Bezier with control points

    # Lagrangian or Bezier with inner points
    if bl.curvetype in (2,3)
        split_curve(bl.coords, bl, msh)
        return
    end

    # Polyline
    if bl.closed && n<3; error("At least thre points are required for closed polyline in BlockInset") end

    bl._endpoint = nothing
    for i=1:n
        if i==n && !bl.closed; break end
        if i<n
            coords = bl.coords[i:i+1,:]
        else
            coords = [ bl.coords[n,:] ; bl.coords[1,:] ]
        end
        split_curve(coords, bl, msh)
    end
end

function get_point(s, coords, curvetype)
    s = s>1.0? 1.0 : s
    if curvetype<=2; 
        return interLagrange(s, coords) 
    else
        return cubicBezier  (s, coords) 
    end
end

function split_curve(coords, bl::BlockInset, msh::Mesh)
    # Constants
    TINY = 1.0e-4
    TOL  = 1.0e-9
    JMP  = 1.0e-2
    nits = 25
    shape   = bl.shape
    npoints = shape==LIN2? 2 : 3
    lnkshape = shape==LIN2? LINK2 : LINK3
    curvetype = bl.curvetype

    # Initial conditions
    bdist = 0.0     # boundary function initial value
    len   = 1.0

    # Defining required vectors
    X1 = vec(coords[  1,:])
    Xn = vec(coords[end,:])

    # Find the initial and final element
    ecells = Array(Cell,0)
    s0 = get_point(TINY,coords,curvetype)
    #ncell::Union(Cell, Nothing); # next cell
    #icell::Union(Cell, Nothing); # initial cell
    icell = find_cell(s0, msh.cells, msh.bins, TOL, ecells) # The first tresspased cell

    if icell == nothing
        error("Inset limits outside the mesh")
    end

    # Initializing more variables
    ccell  = icell
    points = Array(Point, npoints)
    bl._endpoint   = nothing
    end_reached = false
    s  = 0.0
    sp = 0.0
    nits = int(1./JMP)

    # Splitting inset
    k = 0
    while true
        k +=1
        ccell_coords = get_coords(ccell)
        # Default step
        step  = 0.50*(1.0-s)

        # Finding step
        st = s     # trial point
        for i=1:nits
            st += JMP
            if st>1.0; break end
            X = get_point(st, coords, curvetype)
            is_in = is_inside(ccell.shape, ccell_coords, X, TOL)
            if !is_in
                step  = 0.5*(st-s)
                break
            end
        end

        s += step
        X  = get_point(s, coords, curvetype)
        n  = ifloor(log(2, step/TOL)) + 1  # number of required iterations to find intersection

        for i=1:n

            step *= 0.5+TOL
            is_in = is_inside(ccell.shape, ccell_coords, X, TOL)
            if is_in
                s += step
            else
                s -= step
            end

            X = get_point(s, coords, curvetype)

            R     = inverse_map(ccell.shape, ccell_coords, X)
            bdist = bdistance(ccell.shape, R)
        end

        # Check if end was reached
        if s > len - TINY
            end_reached = true
            if curvetype<=1; X=Xn end
            # TODO test (check also incomplete Bezier...)
        end

        # Counter
        bl.icount += end_reached? 0 : 1

        # Getting line cell points
        if bl._endpoint==nothing
            P1 = Point(X1)
            push!(msh.points, P1)
        else
            P1 = bl._endpoint
        end

        if !(bl.closed && end_reached)
            P2 = Point(X) 
            push!(msh.points, P2)
        else
            P2 = bl._startpoint
        end

        if npoints==2
            Ps = [P1, P2]
        else
            P3 = Point(get_point( (sp+s)/2.0, coords, curvetype))
            push!(msh.points, P3)
            Ps = [P1, P2, P3]
        end

        if bl._startpoint == nothing; bl._startpoint = P1 end
        bl._endpoint = P2

        # Saving line cell 
        lcell = Cell(shape, Ps, bl.tag)
        push!(msh.cells, lcell)

        # Create a continuous joint element
        lnkpts  = [ ccell.points, lcell.points ]
        lnkcell = Cell(lnkshape, lnkpts, bl.tag)
        push!(msh.cells, lnkcell)

        ccell.crossed = true

        if end_reached
            return
        end

        # Preparing for the next iteration
        #ecells = [ ccell ]
        ncell  = find_cell(get_point(s + TINY, coords, curvetype), msh.cells, msh.bins, TOL, [ccell])
        if ncell == nothing
            error("Hole found while searching for next tresspassed cell")
        end

        ccell = ncell
        sp = s
        s = s+TINY
    end
end
