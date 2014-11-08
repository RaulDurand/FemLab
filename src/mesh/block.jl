export Block2D, Block3D, BlockTruss
import Base.copy

### Type Block
abstract Block

### Block types:

include("block_inset.jl")


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

copy(bl::BlockTruss) = BlockTruss(copy(bl.coords), bl.conns, shape=bl.shape, tag=bl.tag)



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

copy(bl::Block2D) = Block2D(copy(bl.coords), nx=bl.nx, ny=bl.ny, shape=bl.shape, tag=bl.tag)



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

copy(bl::Block3D) = Block2D(copy(bl.coords), nx=bl.nx, ny=bl.ny, nz=bl.nz, shape=bl.shape, tag=bl.tag)



# Functions for blocks

function move(bl::Block;x=0.0, y=0.0, z=0.0)
    n = size(bl.coords, 1)
    bl.coords[1:n, 1] += x
    bl.coords[1:n, 2] += y
    bl.coords[1:n, 3] += z
    return bl
end


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

                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
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
                    if size(bl.coords,1)==8
                        N = shape_func(HEX8, [r, s, t])
                    else
                        N = shape_func(HEX20, [r, s, t])
                    end
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

                    r = (1.0/nx)*(i-1) - 1.0
                    s = (1.0/ny)*(j-1) - 1.0
                    t = (1.0/nz)*(k-1) - 1.0
                    if size(bl.coords,1)==8
                        N = shape_func(HEX8, [r, s, t])
                    else
                        N = shape_func(HEX20, [r, s, t])
                    end
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

        for k = 1:2:2*nz
            for j = 1:2:2*ny
                for i = 1:2:2*nx


                        #println(" $i  $j  $k ")
                        #aa= p_arr[i  , j  , k  ]
                        #aa= p_arr[i+2, j  , k  ]
                        #aa= p_arr[i+2, j+2, k  ]
                        #aa= p_arr[i  , j+2, k  ]
                        #aa= p_arr[i  , j  , k+2]
                        #aa= p_arr[i+2, j  , k+2]
                        #aa= p_arr[i+2, j+2, k+2]
                        #aa= p_arr[i  , j+2, k+2]
                        #
                        #aa= p_arr[i+1, j  , k  ]
                        #aa= p_arr[i+2, j+1, k  ]
                        #aa= p_arr[i+1, j+2, k  ]
                        #aa= p_arr[i  , j+1, k  ]
                        #aa= p_arr[i+1, j  , k+2]
                        #aa= p_arr[i+2, j+1, k+2]
                        #aa= p_arr[i+1, j+2, k+2]
                        #aa= p_arr[i  , j+1, k+2]
                       #
                        #aa= p_arr[i  , j  , k+1]
                        #aa= p_arr[i+2, j  , k+1]
                        #aa= p_arr[i+2, j+2, k+1]
                        #aa= p_arr[i  , j+2, k+1]
#


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

