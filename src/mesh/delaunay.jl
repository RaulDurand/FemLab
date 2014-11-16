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


type TPoint
    x::Float64
    y::Float64
    id::Int64
    function TPoint(x,y)
        this = new(x,y)
        this.id = -1
        return this
    end
end


type TCell
    points::Array{TPoint,1}
    adjacs::Array{Union(TCell,Nothing),1}

    function TCell(p0, p1, p2, t1=nothing, t2=nothing, t3=nothing)
        this = new()
        this.points = [p0, p1, p2]
        this.adjacs = [t1, t2, t3]
        return this
    end

end


import Base.repr
function repr(cell::TCell)
    repr(cell.points)
end

import Base.show
function show(io::IO, cell::TCell)
    print(repr(cell))
end


inc(v, n) = (v+n-1)%3 + 1

#inc(v, n) = v+n>3? (v+n)%3 : v+n

function find_cell(x, y, cell) # recursive search algorithm
    for i=1:3
        pa = cell.points[i]
        pb = cell.points[inc(i,1)]
        if (pb.x-pa.x)*(y-pa.y) - (pb.y-pa.y)*(x-pa.x) < 0.0 # point lies on right side
            return find_cell(x, y, cell.adjacs[i])
        end
    end
    return cell
end


function swap_cells(cell0, cell1)
    # search for p0 and v0 v1
    local p0, v0, v1
    for i=1:3
        p0 = cell0.points[i]
        if !(p0 in cell1.points) 
            v0 = cell0.points[inc(i,2)]
            v1 = cell0.points[inc(i,1)]
            break 
        end
    end

    #@show cell0.points
    #@show cell1.points

    pos_p0 = findfirst(cell0.points, p0)
    #@show pos_p0

    # Get p1
    #pos_p1 = findfirst(cell1.points, v1)
    pos_p1 = inc(findfirst(cell0.points, v1), 1)
    p1  = cell1.points[pos_p1]

    # Calculate if p0 from cell 0 is in circuncircle of cell1
    x01 = v0.x - p1.x  ;  y01 = v0.y - p1.y
    x11 = v1.x - p1.x  ;  y11 = v1.y - p1.y
    x00 = v0.x - p0.x  ;  y00 = v0.y - p0.y
    x10 = v1.x - p0.x  ;  y10 = v1.y - p0.y

    if (x01*x11 + y01*y11)*(x10*y00 - x00*y10) < (y01*x11 - x01*y11)*(x10*x00 + y00*y10)

        # related adjacent cells
        #@show inc(pos_p1,2)
        a = cell1.adjacs[inc(pos_p1,2)]
        b = cell1.adjacs[pos_p1]
        c = cell0.adjacs[inc(pos_p0,2)]
        d = cell0.adjacs[pos_p0]

        cell0.points = [p0, v1, p1]
        cell1.points = [p0, p1, v0]
        cell0.adjacs = [d, a, cell1]
        cell1.adjacs = [cell0, b, c]

        if a != nothing; a.adjacs[ findfirst(a.points, p1) ] = cell0 end
        if c != nothing; c.adjacs[ findfirst(c.points, p0) ] = cell1 end

        return true
    end

    return false

end



function triangulate(vertices::Array{Float64,2})
    cells  = TCell[]
    points = TPoint[]
    n = size(vertices,1)

    minx, miny = minimum(vertices,1)
    maxx, maxy = maximum(vertices,1)

    Δ = 10.*max(maxx-minx, maxy-miny, 1.0)
    @show Δ 

    pp0 = TPoint(minx-Δ  , miny-Δ)
    pp1 = TPoint(maxx+Δ  , miny-Δ)
    pp2 = TPoint(minx-Δ  , maxy+Δ)
    @show pp0
    @show pp1
    @show pp2

    scell = TCell(pp0, pp1, pp2)
    push!(cells, scell)
    
    for i=1:n
        p = TPoint(vertices[i,1:end]...)
        push!(points, p)

        @show length(cells)
        @show p.x, p.y

        cell = find_cell(p.x, p.y, cells[end])

        # Get cell info
        p0, p1, p2 = cell.points
        a0, a1, a2 = cell.adjacs

        # Define new cells
        t0 = cell
        t0.points = [p, p0, p1]
        t1 = TCell(p, p1, p2); push!(cells, t1)
        t2 = TCell(p, p2, p0); push!(cells, t2)
        t0.adjacs = [t2, a0, t1]
        t1.adjacs = [t0, a1, t2]
        t2.adjacs = [t1, a2, t0]

        # update neighbors in adjacent cells
        if a0 != nothing; a0.adjacs[ findfirst(a0.points, p1) ] = t0 end
        if a1 != nothing; a1.adjacs[ findfirst(a1.points, p2) ] = t1 end
        if a2 != nothing; a2.adjacs[ findfirst(a2.points, p0) ] = t2 end

        toswap = [t0, t1, t2]

        while length(toswap)>0
            cell0 = pop!(toswap)
            cell1 = cell0.adjacs[2] # adjacent cell opposito to point P in cell0 cell
            if cell1==nothing
                continue
            end

            if swap_cells(cell0, cell1)
                push!(toswap, cell0)
                push!(toswap, cell1)
            end

        end

    end


    # remove extra cells related to supertriangle
    spoints =  [pp0, pp1, pp2]
    n = length(cells)
    icells = Int64[]
    @show spoints
    for i=1:n
        cell = cells[i] # save positions
        # remove adjacents in remaining cells
        @show cell
        if (pp0 in cell.points) || (pp1 in cell.points) || (pp2 in cell.points)
            for (p,adj) in zip(cell.points, cell.adjacs)
                if adj!=nothing
                    adj.adjacs[ inc(findfirst(adj.points, p), 2) ] = nothing
                end
            end
        else
            push!(icells, i)
        end
    end

    @show icells
    cells = cells[icells]

    # numbering points
    for (i,point) in enumerate(points)
        point.id = i
    end

    # connectivities
    n = length(cells)
    C = Array(Int64, n, 3)
    for (i,cell) in enumerate(cells)
        C[i,1] = cell.points[1].id
        C[i,2] = cell.points[2].id
        C[i,3] = cell.points[3].id
    end

    @show points
    edges = Set{(Int64, Int64)}()
    for (i,cell) in enumerate(cells)
        for (pos1, pos2) in ((1,2), (2,3), (3,1))
            id1, id2 = cell.points[pos1].id, cell.points[pos2].id
            @show id1, id2
            push!(edges, (min(id1,id2), max(id1, id2)))
        end
    end

    n = length(edges)
    E = Array(Int64, n, 2)
    #@show points
    @show edges
    for (i,edge) in enumerate(edges)
        E[i,1] = edge[1]
        E[i,2] = edge[2]
    end

    @show E

    return C, E

end


V = [ 1. 2; 2 1; 1 3 ; 2 4; 3 2; 4 3 ]
C, E = triangulate(V)


using FemLab

bl = BlockTruss(V, E)
mesh = generate_mesh(bl)
save(mesh, "tri.vtk")




