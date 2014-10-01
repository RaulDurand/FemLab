
#using Shape
export Point, Cell, hash, get_coords, get_point, filter, get_faces

### Type Point definition

type Point
    x::Float64
    y::Float64
    z::Float64
    tag::String
    id::Int64
    function Point(x,y,z; tag="")
        const NDIG = 8
        x = round(x, NDIG)
        y = round(y, NDIG)
        z = round(z, NDIG)
        return new(x,y,z,tag,-1)
    end
    function Point(C::Array{Float64,1})
        return Point(C[1], C[2], C[3])
    end
end


### Point methods

import Base.hash
#hash(p::Point) = round(p.x*1001 + p.y*10000001 + p.z*100000000001)
hash(p::Point) = hash((p.x, p.y, p.z))

getcoords(p::Point) = [p.x, p.y, p.z]

function get_point(points::Dict{Uint64,Point}, C::Array{Float64,1})
    hs = hash(Point(C))
    get(points, hs, nothing)
end

function filter(points::Array{Point,1}; x=NaN, y=NaN, z=NaN)
    result = Array(Point, 0)

    for p in points
        if (x==p.x || isnan(x)) && (y==p.y || isnan(y)) && (z==p.z || isnan(z)) 
            push!(result, p)
        end
    end
    return result

end


### Type Cell definition


type Cell
    shape ::ShapeType
    points::Array{Point,1}
    tag ::String
    id  ::Integer
    ndim::Integer
    crossed::Bool
    ocell::Union(Cell,Nothing)
    function Cell(shape::ShapeType, points, tag="", ocell=nothing)
        this = new(shape, points, tag, -1)
        this.ndim = 0
        this.crossed = false
        this.ocell   = ocell
        return this
    end
end


### Cell methods


hash(c::Cell) = sum([ hash(p) for p in c.points])

function get_coords(c::Cell)
    n = length(c.points)
    C = Array(Float64, n, 3)
    for (i,p) in enumerate(c.points)
        C[i,1] = p.x
        C[i,2] = p.y
        C[i,3] = p.z
    end
    return C
end

type Bins
    bins::Array{Array{Cell,1},3}
    bbox::Array{Float64,2}
    lbin::Float64
    function Bins()
        this = new()
        this.bins = Array(Array{Cell,1}, 0, 0, 0)
        this.bbox = zeros(0,0)
        return this
    end
end

function bounding_box(cells::Array{Cell,1})
    # Get all points
    pointsd = Dict{Uint64, Point}()
    for cell in cells
        for point in cell.points
            pointsd[hash(point)] = point
        end
    end
    points = values(pointsd)

    # Get bounding box
    minx = miny = minz =  Inf
    maxx = maxy = maxz = -Inf
    for point in points
        if point.x<minx; minx = point.x end
        if point.y<miny; miny = point.y end
        if point.z<minz; minz = point.z end
        if point.x>maxx; maxx = point.x end
        if point.y>maxy; maxy = point.y end
        if point.z>maxz; maxz = point.z end
    end
    return [ minx miny minz; maxx maxy maxz ]
end

function bounding_box(cell::Cell)
    minx = miny = minz =  Inf
    maxx = maxy = maxz = -Inf
    for point in cell.points
        if point.x<minx; minx = point.x end
        if point.y<miny; miny = point.y end
        if point.z<minz; minz = point.z end
        if point.x>maxx; maxx = point.x end
        if point.y>maxy; maxy = point.y end
        if point.z>maxz; maxz = point.z end
    end
    return [ minx miny minz; maxx maxy maxz ]
end

function build_bins(cells::Array{Cell,1}, bins::Bins)
    # Get all points
    bins.bbox = bounding_box(cells)
    minx, miny, minz = bins.bbox[1,:]

    # Get max global lengths
    Lx, Ly, Lz = bins.bbox[2,:] - bins.bbox[1,:]
    max_L = max(Lx, Ly, Lz)

    # Get max cell lengths
    max_l = 0.0
    for cell in cells
        if !is_solid(cell.shape); continue end
        bbox = bounding_box(cell)
        max  = maximum(bbox[2,:] - bbox[1,:])
        #if max>0.9
            #@out bbox
            #@out max
            #@out cell.id
            #@out get_coords(cell)
            #exit()
        #end
        if max>max_l; max_l = max end
    end

    # Get number of divisions
    ndiv = min(50, 1*ifloor(max_L/max_l)) # calibrate for bins efficiency

    lbin = max_L/ndiv     # Get bin length
    bins.lbin = lbin

    nx = ifloor(Lx/lbin) + 1
    ny = ifloor(Ly/lbin) + 1
    nz = ifloor(Lz/lbin) + 1

    #@out max_L
    #@out max_l
    #@out ndiv
    #@out (nx,ny,nz)
    #exit()

    # Allocate bins
    bins.bins = Array(Array{Cell,1}, nx, ny, nz)
    for k=1:nz, j=1:ny, i=1:nx
        bins.bins[i,j,k] = Array(Cell,0)
    end

    # Fill bins
    for cell in cells
        bbox = bounding_box(cell)
        X, Y, Z  = bbox[:,1], bbox[:,2], bbox[:,3]
        verts    = [ [x y z] for z in Z, y in Y, x in X ]
        cell_pos = Set()
        for V in verts
            x, y, z = V
            ix = ifloor((x - minx)/lbin) + 1
            iy = ifloor((y - miny)/lbin) + 1
            iz = ifloor((z - minz)/lbin) + 1
            push!(cell_pos, (ix, iy, iz))
        end

        for (ix, iy, iz) in cell_pos
            push!(bins.bins[ix, iy, iz], cell)
        end
    end
end

function find_cell(X::Array{Float64,1}, cells::Array{Cell,1}, bins::Bins, Tol::Float64, exc_cells::Array{Cell,1})
    # Point coordinates
    x, y, z = [X, 0][1:3]
    lbin = bins.lbin

    # Build bins if empty
    if length(bins.bins) == 0
        build_bins(cells, bins)
    end

    for attempt=1:2 
        Cmin = reshape(bins.bbox[1,:],3)
        Cmax = reshape(bins.bbox[2,:],3)
        lbin = bins.lbin

        if any(X.<Cmin-Tol) || any(X.>Cmax+Tol)
            error("point outside bounding box")
        end

        # Find bin index
        ix = ifloor((x - Cmin[1])/lbin) + 1
        iy = ifloor((y - Cmin[2])/lbin) + 1
        iz = ifloor((z - Cmin[3])/lbin) + 1

        # Search cell in bin
        bin = bins.bins[ix, iy, iz]
        for cell in bin
            coords = get_coords(cell)
            if is_inside(cell.shape, coords, X, Tol) && !(cell in exc_cells)
                return cell
            end
        end

        # If not found in the first try then rebuild bins
        if attempt==1
            println("Bin search failed. Rebuilding bins...")
            build_bins(cells, bins) 
        end
    end

    println("Bin search failed")
    return nothing
end

FACETS_IDXS = [
    TRI3   => ( (1, 2),                   (2, 3),                    (3, 1)                                                                                                  ),
    TRI6   => ( (1, 2, 4),                (2, 3, 5),                 (3, 1, 6)                                                                                               ),
    TRI9   => ( (1, 2, 4, 7),             (2, 3, 5, 8),              (3, 1, 6, 9)                                                                                            ),
    QUAD4  => ( (1, 2),                   (2, 3),                    (3, 4),                   (4, 1)                                                                        ),
    QUAD8  => ( (1, 2, 5),                (2, 3, 6),                 (3, 4, 7),                (4, 1, 8)                                                                     ),
    QUAD16 => ( (1, 2, 5),                (2, 3, 6),                 (3, 4, 7),                (4, 1, 8)                                                                     ),
    TET4   => ( (1, 4, 3),                (1, 2, 4),                 (1, 3, 2),                (2, 3, 4)                                                                     ),
    TET10  => ( (1, 4, 3, 8, 10, 7),      (1, 2, 4, 5, 9, 8),        (1, 3, 2, 7, 6, 5),       (2, 3, 4, 6, 10, 9)                                                           ),
    HEX8   => ( (1, 5, 8, 4),             (2, 3, 7, 6),              (1, 2, 6, 5),             (3, 4, 8, 7),             (1, 4, 3, 2),              (5, 6, 7, 8)             ),
    HEX20  => ( (1, 5, 8, 4,17,16,20,12), (2, 3, 7, 6, 10,19,14,18), (1, 2, 6, 5, 9,18,13,17), (3, 4, 8, 7,11,20,15,19), (1, 4, 3, 2,12,11, 10, 9), (5, 6, 7, 8,13,14,15,16) ),
    ]

FACETS_SHAPE = [
    QUAD4  => LIN2 ,
    TRI3   => LIN2 ,
    TRI6   => LIN3 ,
    TRI9   => LIN4 ,
    QUAD4  => LIN2 ,
    QUAD8  => LIN3 ,
    QUAD16 => LIN4 ,
    TET4   => TRI3 ,
    TET10  => TRI6 ,
    HEX8   => QUAD4,
    HEX20  => QUAD8
]

function get_faces(cell::Cell)
    faces  = Array(Cell,0)
    if !haskey(FACETS_IDXS, cell.shape)
        return faces
    end

    f_idxs = FACETS_IDXS[cell.shape]

    for f_idx in f_idxs
        points = [ cell.points[i] for i in f_idx]
        face   = Cell(FACETS_SHAPE[cell.shape], points, cell.tag, cell)
        push!(faces, face)
    end

    return faces
end

#end #module
