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

export BlockInset

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
        
        if size(coords,2)==2
            nrows = size(coords,1)
            coords = [ coords  zeros(nrows)]
        end

        this = new(coords, ctype, closed, shape, tag, id)
        this.icount = 0
        this._endpoint   = nothing
        this._startpoint = nothing
        return this
    end
end

copy(bl::BlockInset) = BlockInset(copy(bl.coords), curvetype=bl.curvetype, closed=bl.closed, shape=bl.shape, tag=bl.tag)


function cubicBezier(s::Float64, PorQ::Array{Float64,2}, isQ::Bool=false)
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

function interLagrange(s::Float64, coords::Array{Float64,2})
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
        split_curve(bl.coords, bl, false, msh)
        return
    end

    # Polyline
    if bl.closed && n<3; error("At least three points are required for closed polyline in BlockInset") end

    bl._endpoint = nothing
    for i=1:n-1
        coords = bl.coords[i:i+1,:]
        split_curve(coords, bl, false, msh)
    end
    if bl.closed
        coords = [ bl.coords[n,:] ; bl.coords[1,:] ]
        split_curve(coords, bl, true, msh)
    end
end

function get_point(s::Float64, coords::Array{Float64,2}, curvetype::ShapeType)
    s = s>1.0? 1.0 : s
    if curvetype<=2; 
        return interLagrange(s, coords) 
    else
        return cubicBezier  (s, coords) 
    end
end

function split_curve(coords::Array{Float64,2}, bl::BlockInset, closed::Bool, msh::Mesh)
    # Constants
    TINY = 1.0e-4
    TOL  = 1.0e-9
    JMP  = 1.0e-2


    #TINY = 1.0E-3 ##
    #TOL  = 1.0E-4 ##
    #JMP  = 1.0 ##

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
    icell = find_cell(s0, msh.cells, msh.bins, TOL, ecells) # The first tresspased cell

    if icell == nothing
        error("Inset limits outside the mesh")
    end

    # Initializing more variables
    ccell  = icell
    points = Array(Point, npoints)
    bl._endpoint = nothing
    end_reached  = false
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

        # R     = inverse_map(ccell.shape, ccell_coords, X) ##
        # bdist = bdistance(ccell.shape, R) ##
        #@printf(" %d  &  %10.6f  &  %10.6f  & %10.6f  & %10.6f \\\\\n", 0, s, X[1], X[2], bdist)

        for i=1:n

            #step *= 0.5+TOL
            step *= 0.5
            is_in = is_inside(ccell.shape, ccell_coords, X, TOL)
            if is_in
                s += step
            else
                s -= step
            end

            X = get_point(s, coords, curvetype)

            R     = inverse_map(ccell.shape, ccell_coords, X)
            bdist = bdistance(ccell.shape, R)
            #@printf(" %d  &  %10.6f  &  %10.6f  & %10.6f  & %10.6f \\\\\n", i, s, X[1], X[2], bdist)
        end


        # println()
        # ds = (1-s)/2
        # @printf(" %d  &  %10.6f  &  %2d  &  %10.6f  & %10.6f  & %10.6f \\\\\n", 1, ds, n, s, X[1], X[2])
        # println()

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

        if !(closed && end_reached)
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
