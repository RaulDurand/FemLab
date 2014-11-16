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


include("quadrature.jl")

export get_vtk_type, get_shape_from_vtk
export LIN2, LIN3, TRI3, TRI6, QUAD4, QUAD8, TET4, TET10, HEX8, HEX20 
export LINK1, LINK2, LINK3, LIN4, TRI9, TRI10, QUAD9, QUAD12, QUAD16
export ShapeType
export shape_func, deriv_func, local_coords
export get_ip_coords
export is_solid, is_line, is_joint, is_inside
export inverse_map, extrapolator

# Shape types:
const LIN2  = 3
const LIN3  = 21
const TRI3  = 5
const TRI6  = 22
const QUAD4 = 9
const QUAD8 = 23
const TET4  = 10
const TET10 = 24
const HEX8  = 12
const HEX20 = 25

const LIN4   = 31
const TRI9   = 33
const TRI10  = 34
const QUAD9  = 37
const QUAD12 = 38
const QUAD16 = 39

const LINK1  = 51
const LINK2  = 52
const LINK3  = 53

typealias ShapeType Integer

function get_vtk_type(shape::ShapeType)
    if shape in [LINK1, LINK2, LINK3, LIN4, TRI9, TRI10, QUAD9, QUAD12, QUAD16]
        return 2 # vtk_poly_vertex
    end

    return shape # Conventional:
end


function get_shape_from_vtk(vtk_type::Int64, npoints=None)
    if vtk_type!=2 return vtk_type end # poly_vertex

    if     npoints==1   LINK1
    elseif npoints==2   LINK2
    elseif npoints==3   LINK3
    elseif npoints==9   TRI9
    elseif npoints==10  TRI10
    elseif npoints==12  QUAD12
    elseif npoints==16  QUAD16
    end
end


coords_lin2 =
[ -1.0  1.0 
   1.0  1.0 ]

coords_tri6 =
[ 0.0  0.0  1.0
  1.0  0.0  1.0
  0.0  1.0  1.0 
  0.5  0.0  1.0 
  0.5  0.5  1.0 
  0.0  0.5  1.0 ]


_1_3 = 1.0/3.0; _2_3 = 2.0/3.0

coords_tri9 =
[ 0.0   0.0  1.0
  1.0   0.0  1.0
  0.0   1.0  1.0
  
  _1_3   0.0  1.0
  _2_3  _1_3  1.0
   0.0  _2_3  1.0
  
  _2_3   0.0  1.0
  _1_3  _2_3  1.0
   0.0  _1_3  1.0 ]

coords_tri10 =
[ 0.0   0.0  1.0 
  1.0   0.0  1.0 
  0.0   1.0  1.0 
  
  _1_3   0.0  1.0 
  _2_3  _1_3  1.0 
   0.0  _2_3  1.0 
  
  _2_3   0.0  1.0 
  _1_3  _2_3  1.0 
   0.0  _1_3  1.0 
  
  _1_3  _1_3  1.0 ]


coords_quad4 = 
[ -1.0 -1.0  1.0
   1.0 -1.0  1.0
   1.0  1.0  1.0
  -1.0  1.0  1.0 ]


immutable type Typed{N} end

shape_func(shape::ShapeType, R::Array{Float64,1}) = shape_func(Typed{shape}(), R)
deriv_func(shape::ShapeType, R::Array{Float64,1}) = deriv_func(Typed{shape}(), R)

function shape_func(::Typed{LIN2}, R::Array{Float64,1})
    r = R[1]
    N = Array(Float64, 2)
    N[1] = 0.5*(1-r)
    N[2] = 0.5*(1+r)
    return N
end

function deriv_func(::Typed{LIN2}, R::Array{Float64,1})
    D = Array(Float64, 1, 2)
    D[1,1] = -0.5
    D[1,2] =  0.5
    return D
end

function shape_func(::Typed{LIN3}, R::Array{Float64,1})
    r = R[1]
    N = Array(Float64, 3)
    N[1] = 0.5*(r*r - r)
    N[2] = 0.5*(r*r + r)
    N[3] = 1.0 - r*r
    return N
end

function deriv_func(::Typed{LIN3}, R::Array{Float64,1})
    r = R[1]
    D = Array(Float64, 1, 3)
    D[1, 1] = r - 0.5
    D[1, 2] = r + 0.5
    D[1, 3] = -2.0*r
    return D
end

function shape_func(::Typed{QUAD4}, R::Array{Float64,1})
    #     4                        3
    #       @--------------------@
    #       |               (1,1)|
    #       |       s ^          |
    #       |         |          |
    #       |         |          |
    #       |         +----> r   |
    #       |       (0,0)        |
    #       |                    |
    #       |                    |
    #       |(-1,-1)             |
    #       @--------------------@
    #     1                        2
    #
    r, s = R[1:2]
    N = Array(Float64,4)
    N[1] = 0.25*(1.0-r-s+r*s)
    N[2] = 0.25*(1.0+r-s-r*s)
    N[3] = 0.25*(1.0+r+s+r*s)
    N[4] = 0.25*(1.0-r+s-r*s)
    return N
end

function deriv_func(::Typed{QUAD4}, R::Array{Float64,1})
    r, s = R[1:2]
    D = Array(Float64, 2, 4)
    D[1,1] = 0.25*(-1.0+s);   D[2,1] = 0.25*(-1.0+r)
    D[1,2] = 0.25*(+1.0-s);   D[2,2] = 0.25*(-1.0-r)
    D[1,3] = 0.25*(+1.0+s);   D[2,3] = 0.25*(+1.0+r)
    D[1,4] = 0.25*(-1.0-s);   D[2,4] = 0.25*(+1.0-r)
    return D
end

function shape_func(::Typed{QUAD8}, R::Array{Float64,1})
    #     4           7            3
    #       @---------@----------@
    #       |               (1,1)|
    #       |       s ^          |
    #       |         |          |
    #       |         |          |
    #     8 @         +----> r   @ 6
    #       |       (0,0)        |
    #       |                    |
    #       |                    |
    #       |(-1,-1)             |
    #       @---------@----------@
    #     1           5            2
    #
    r, s = R[1:2]
    N = Array(Float64,8)
    rp1=1.0+r; rm1=1.0-r;
    sp1=1.0+s; sm1=1.0-s;
    N[1] = 0.25*rm1*sm1*(rm1+sm1-3.0)
    N[2] = 0.25*rp1*sm1*(rp1+sm1-3.0)
    N[3] = 0.25*rp1*sp1*(rp1+sp1-3.0)
    N[4] = 0.25*rm1*sp1*(rm1+sp1-3.0)
    N[5] = 0.50*sm1*(1.0-r*r)
    N[6] = 0.50*rp1*(1.0-s*s)
    N[7] = 0.50*sp1*(1.0-r*r)
    N[8] = 0.50*rm1*(1.0-s*s)
    return N
end

function deriv_func(::Typed{QUAD8}, R::Array{Float64,1})
    r, s = R[1:2]
    D = Array(Float64, 2, 8)
    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s

    D[1,1] = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0)
    D[1,2] =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0)
    D[1,3] =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0)
    D[1,4] = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0)
    D[1,5] = - r * sm1
    D[1,6] =   0.50 * (1.0 - s * s)
    D[1,7] = - r * sp1
    D[1,8] = - 0.5 * (1.0 - s * s)

    D[2,1] = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0)
    D[2,2] = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0)
    D[2,3] =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0)
    D[2,4] =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0)
    D[2,5] = - 0.50 * (1.0 - r * r)
    D[2,6] = - s * rp1
    D[2,7] =   0.50 * (1.0 - r * r)
    D[2,8] = - s * rm1
    return D
end

function shape_func(::Typed{HEX8}, R::Array{Float64,1})
    # Local IDs
    #                  Nodes                                   Faces
    #     z
    #     |           5                  8
    #    ,+--y         @________________@                    +________________+
    #  x'            ,'|              ,'|                  ,'|              ,'|
    #              ,'  |            ,'  |                ,'  |  ___       ,'  |
    #            ,'    |          ,'    |              ,'    |,'5,'  [0],'    |
    #      6   ,'      |      7 ,'      |            ,'      |~~~     ,'      |
    #        @'===============@'        |          +'===============+'  ,'|   |
    #        |         |      |         |          |   ,'|   |      |   |3|   |
    #        |         |      |         |          |   |2|   |      |   |,'   |
    #        |       1 @______|_________@          |   |,'   +______|_________+
    #        |       ,'       |       ,' 4         |       ,'       |       ,'
    #        |     ,'         |     ,'             |     ,' [1]  ___|     ,'
    #        |   ,'           |   ,'               |   ,'      ,'4,'|   ,'
    #        | ,'             | ,'                 | ,'        ~~~  | ,'
    #        @________________@'                   +________________+'
    #      2                   3

    r, s, t = R[1:3]
    N = Array(Float64,8)
    N[1] = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t)
    N[2] = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t)
    N[3] = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t)
    N[4] = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t)
    N[5] = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t)
    N[6] = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t)
    N[7] = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t)
    N[8] = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t)
    return N
end

function deriv_func(::Typed{HEX8}, R::Array{Float64,1})
    r, s, t = R
    st = s*t
    rt = r*t
    rs = r*s
    D = Array(Float64, 3, 8)
    D[1,1] = -1.0+s+t-st;   D[2,1]=-1.0+r+t-rt;   D[3,1]=-1.0+r+s-rs
    D[1,2] = +1.0-s-t+st;   D[2,2]=-1.0-r+t+rt;   D[3,2]=-1.0-r+s+rs
    D[1,3] = +1.0+s-t-st;   D[2,3]=+1.0+r-t-rt;   D[3,3]=-1.0-r-s-rs
    D[1,4] = -1.0-s+t+st;   D[2,4]=+1.0-r-t+rt;   D[3,4]=-1.0+r-s+rs
    D[1,5] = -1.0+s-t+st;   D[2,5]=-1.0+r-t+rt;   D[3,5]=+1.0-r-s+rs
    D[1,6] = +1.0-s+t-st;   D[2,6]=-1.0-r-t-rt;   D[3,6]=+1.0+r-s-rs
    D[1,7] = +1.0+s+t+st;   D[2,7]=+1.0+r+t+rt;   D[3,7]=+1.0+r+s+rs
    D[1,8] = -1.0-s-t-st;   D[2,8]=+1.0-r+t-rt;   D[3,8]=+1.0-r+s-rs
    D = 0.125*D

    return D
end

function shape_func(::Typed{HEX20}, R::Array{Float64,1})
    # Local IDs
    #                   Vertices                               Faces
    #     t
    #     |           5        16        8
    #    ,+--s         @-------@--------@                   +----------------+
    #  r'            ,'|              ,'|                 ,'|              ,'|
    #           13 @'  |         15 ,'  |               ,'  |  ___       ,'  |
    #            ,'    |17        ,@    |20           ,'    |,'6,'  [1],'    |
    #      6   ,'      @      7 ,'      @           ,'      |~~~     ,'      |
    #        @'=======@=======@'        |         +'===============+'  ,'|   |
    #        |      14 |      |         |         |   ,'|   |      |   |4|   |
    #        |         |      |  12     |         |   |3|   |      |   |,'   |
    #     18 |       1 @- - - | @- - - -@         |   |,'   +- - - | +- - - -+
    #        @       ,'       @       ,' 4        |       ,'       |       ,'
    #        |   9 @'      19 |     ,'            |     ,' [2]  ___|     ,'
    #        |   ,'           |   ,@ 11           |   ,'      ,'5,'|   ,'
    #        | ,'             | ,'                | ,'        ~~~  | ,'
    #        @-------@--------@'                  +----------------+'
    #      2        10         3

    r, s, t = R

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    N = Array(Float64,20)
    N[ 1] = 0.125*rm1*sm1*tm1*(-r-s-t-2.0)
    N[ 2] = 0.125*rp1*sm1*tm1*( r-s-t-2.0)
    N[ 3] = 0.125*rp1*sp1*tm1*( r+s-t-2.0)
    N[ 4] = 0.125*rm1*sp1*tm1*(-r+s-t-2.0)
    N[ 5] = 0.125*rm1*sm1*tp1*(-r-s+t-2.0)
    N[ 6] = 0.125*rp1*sm1*tp1*( r-s+t-2.0)
    N[ 7] = 0.125*rp1*sp1*tp1*( r+s+t-2.0)
    N[ 8] = 0.125*rm1*sp1*tp1*(-r+s+t-2.0)
    N[ 9] = 0.25*(1.0-r*r)*sm1*tm1
    N[10] = 0.25*rp1*(1.0-s*s)*tm1
    N[11] = 0.25*(1.0-r*r)*sp1*tm1
    N[12] = 0.25*rm1*(1.0-s*s)*tm1
    N[13] = 0.25*(1.0-r*r)*sm1*tp1
    N[14] = 0.25*rp1*(1.0-s*s)*tp1
    N[15] = 0.25*(1.0-r*r)*sp1*tp1
    N[16] = 0.25*rm1*(1.0-s*s)*tp1
    N[17] = 0.25*rm1*sm1*(1.0-t*t)
    N[18] = 0.25*rp1*sm1*(1.0-t*t)
    N[19] = 0.25*rp1*sp1*(1.0-t*t)
    N[20] = 0.25*rm1*sp1*(1.0-t*t)
    return N
end

function deriv_func(::Typed{HEX20}, R::Array{Float64,1})
    r, s, t = R

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    D = Array(Float64, 3, 20)
    # Derivatives with respect to r
    D[1, 1] = -0.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[1, 2] =  0.125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1
    D[1, 3] =  0.125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[1, 4] = -0.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[1, 5] = -0.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[1, 6] =  0.125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[1, 7] =  0.125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[1, 8] = -0.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1
    D[1, 9] = -0.5*r*sm1*tm1
    D[1,10] =  0.25*(1-s*s)*tm1
    D[1,11] = -0.5*r*sp1*tm1
    D[1,12] = -0.25*(1-s*s)*tm1
    D[1,13] = -0.5*r*sm1*tp1
    D[1,14] =  0.25*(1-s*s)*tp1
    D[1,15] = -0.5*r*sp1  *tp1
    D[1,16] = -0.25*(1-s*s)*tp1
    D[1,17] = -0.25*sm1*(1-t*t)
    D[1,18] =  0.25*sm1*(1-t*t)
    D[1,19] =  0.25*sp1*(1-t*t)
    D[1,20] = -0.25*sp1*(1-t*t)

    # Derivatives with respect to s
    D[2, 1] = -0.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[2, 2] = -0.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[2, 3] =  0.125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[2, 4] =  0.125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1
    D[2, 5] = -0.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[2, 6] = -0.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1
    D[2, 7] =  0.125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[2, 8] =  0.125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[2, 9] = -0.25*(1-r*r)*tm1
    D[2,10] = -0.5*s*rp1*tm1
    D[2,11] =  0.25*(1-r*r)*tm1
    D[2,12] = -0.5*s*rm1*tm1
    D[2,13] = -0.25*(1-r*r)*tp1
    D[2,14] = -0.5*s*rp1*tp1
    D[2,15] =  0.25*(1-r*r)*tp1
    D[2,16] = -0.5*s*rm1*tp1
    D[2,17] = -0.25*rm1*(1-t*t)
    D[2,18] = -0.25*rp1*(1-t*t)
    D[2,19] =  0.25*rp1*(1-t*t)
    D[2,20] =  0.25*rm1*(1-t*t)

    # Derivatives with respect to t
    D[3, 1] = -0.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[3, 2] = -0.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[3, 3] = -0.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1
    D[3, 4] = -0.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[3, 5] =  0.125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1
    D[3, 6] =  0.125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[3, 7] =  0.125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[3, 8] =  0.125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[3, 9] = -0.25*(1-r*r)*sm1
    D[3,10] = -0.25*rp1*(1-s*s)
    D[3,11] = -0.25*(1-r*r)*sp1
    D[3,12] = -0.25*rm1*(1-s*s)
    D[3,13] =  0.25*(1-r*r)*sm1
    D[3,14] =  0.25*rp1*(1-s*s)
    D[3,15] =  0.25*(1-r*r)*sp1
    D[3,16] =  0.25*rm1*(1-s*s)
    D[3,17] = -0.5*t*rm1*sm1
    D[3,18] = -0.5*t*rp1*sm1
    D[3,19] = -0.5*t*rp1*sp1
    D[3,20] = -0.5*t*rm1*sp1

    return D
end


function local_coords(st::ShapeType)
    if     sh == LIN2   coords_lin2
    elseif sh == QUAD4  coords_quad4
    end
end


IP_FEM = [
    LIN2    => [0 => LIN_IP2,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ],
    LIN3    => [0 => LIN_IP3,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ],
    LIN4    => [0 => LIN_IP3,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ],
    TRI3    => [0 => TRI_IP1,  1 => TRI_IP1,  3 => TRI_IP3,   6 => TRI_IP6                   ],
    TRI6    => [0 => TRI_IP3,  3 => TRI_IP3,  6 => TRI_IP6                                   ],
    TRI9    => [0 => TRI_IP6,  3 => TRI_IP3,  6 => TRI_IP6                                   ],
    TRI10   => [0 => TRI_IP6,  3 => TRI_IP3,  6 => TRI_IP6                                   ],
    LINK2   => [0 => LIN_IP2,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ],
    LINK3   => [0 => LIN_IP3,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ],
    QUAD4   => [0 => QUAD_IP2, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ],
    QUAD8   => [0 => QUAD_IP3, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ],
    QUAD12  => [0 => QUAD_IP4, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ],
    QUAD16  => [0 => QUAD_IP4, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ],
    TET4    => [0 => TET_IP4,  1 => TET_IP1,  4 => TET_IP4,   5 => TET_IP5,  11 => TET_IP11  ],
    TET10   => [0 => TET_IP4,  1 => TET_IP1,  4 => TET_IP4,   5 => TET_IP5,  11 => TET_IP11  ],
    HEX8    => [0 => HEX_IP2,  8 => HEX_IP2, 27 => HEX_IP3                                   ],
    HEX20   => [0 => HEX_IP3,  8 => HEX_IP2, 27 => HEX_IP3                                   ]
]

function get_ip_coords(shape::ShapeType, nips=0)
    if !haskey(IP_FEM, shape)
        error("Ip coordinates for shape ($shape) is not available")
    end
    all_shape_coord = IP_FEM[shape]
    if !haskey(all_shape_coord, nips)
        error("Number of ips ($nips) for shape ($shape) is not available")
    end
    all_shape_coord[nips]
end

function bdistance(shape, R)
    # Returns a real value which is a pseudo distance from a point to the border of an element
    # Arguments:
    #     R - a vector containing the point coordinates
    # Returns:
    #     a real value: if possitive then the point is inside the element and negative otherwise
    
    r, s, t = R
    if shape == TRI3 ;  return min(r, s, 1.0-r-s) end
    if shape == TRI6 ;  return min(r, s, 1.0-r-s) end
    if shape == TRI9 ;  return min(r, s, 1.0-r-s) end
    if shape == TRI10;  return min(r, s, 1.0-r-s) end
    if shape == QUAD4;  return min(1.0 - abs(r), 1.0 - abs(s)) end
    if shape == QUAD8;  return min(1.0 - abs(r), 1.0 - abs(s)) end
    if shape == QUAD12; return min(1.0 - abs(r), 1.0 - abs(s)) end
    if shape == QUAD16; return min(1.0 - abs(r), 1.0 - abs(s)) end
    if shape == TET4 ;  return min(r, s, t, 1.0-r-s-t) end
    if shape == TET10;  return min(r, s, t, 1.0-r-s-t) end
    if shape == HEX8 ;  return min(1.0 - abs(r), 1.0 - abs(s), 1.0 - abs(t)) end
    if shape == HEX20;  return min(1.0 - abs(r), 1.0 - abs(s), 1.0 - abs(t)) end
    error("No boundary distance for shape ($shape)")
end

function get_ndim(shape)
    # Returns the local dimension based on the shape geometry.
    # It does not match necessarily the space where the shape is used.
    
    if shape == LIN2  ; return 1 end
    if shape == LIN3  ; return 1 end
    if shape == LIN4  ; return 1 end
    if shape == LINK2 ; return 1 end
    if shape == LINK3 ; return 1 end
    if shape == TRI3  ; return 2 end
    if shape == TRI6  ; return 2 end
    if shape == TRI9  ; return 2 end
    if shape == TRI10 ; return 2 end
    if shape == QUAD4 ; return 2 end
    if shape == QUAD8 ; return 2 end
    if shape == QUAD12; return 2 end
    if shape == QUAD16; return 2 end
    if shape == TET4  ; return 3 end
    if shape == TET10 ; return 3 end
    if shape == HEX8  ; return 3 end
    if shape == HEX20 ; return 3 end
    error("Unknown shape ($shape)")
end

function get_nnodes(shape)
    # Returns the local dimension based on the shape geometry.
    # It does not match necessarily the space where the shape is used.
    
    if shape == LIN2  ; return  2 end
    if shape == LIN3  ; return  3 end
    if shape == LIN4  ; return  4 end
    if shape == LINK2 ; return  2 end
    if shape == LINK3 ; return  3 end
    if shape == TRI3  ; return  3 end
    if shape == TRI6  ; return  6 end
    if shape == TRI9  ; return  9 end
    if shape == TRI10 ; return 10 end
    if shape == QUAD4 ; return  4 end
    if shape == QUAD8 ; return  8 end
    if shape == QUAD12; return 12 end
    if shape == QUAD16; return 16 end
    if shape == TET4  ; return  4 end
    if shape == TET10 ; return 10 end
    if shape == HEX8  ; return  8 end
    if shape == HEX20 ; return 20 end
    error("Unknown shape ($shape)")
end


function is_line(shape::ShapeType)
    shape in (LIN2, LIN3, LIN4) 
end

function is_joint(shape::ShapeType)
    shape>50? true : false
end

function is_solid(shape::ShapeType)
    # exclude line and link (contact) cells
    !is_line(shape) && !is_joint(shape) ? true : false
end


function inverse_map(shape, coords, X0, Tol=1.0e-7)
    MAXIT = 25
    ndim  = get_ndim(shape)
    R = zeros(ndim)
    C = size(coords,1)==2? [coords zeros(size(coords,1))] : coords
    X = size(X0,1)==2? [X0, 0.0] : X0

    k = 0
    for k=1:MAXIT
        # calculate Jacobian
        D = deriv_func(shape, R)
        J = D*C

        # calculate trial of real coordinates
        N  = shape_func(shape, R)
        Xt = C'*N # interpolating

        # calculate the error
        ΔX = Xt - X
        ΔR = pinv(J)'*ΔX

        # updating local coords R
        R -= ΔR
        if norm(ΔX) < Tol; break end
    end

    #@out k
    if k==MAXIT; print("Warning: max iterations reached in inverse mapping") end

    if ndim==2
        return [ R, 0.0 ]
    else
        return R
    end
end


function is_inside(shape::ShapeType, C::Array{Float64,2}, X::Array{Float64,1}, Tol = 1.e-7)
    if !is_solid(shape) return false end

    # Testing with bounding box
    ndim = size(C,1)
    Cmin = vec(minimum(C,1))
    Cmax = vec(maximum(C,1))
    maxl = maximum(Cmax-Cmin)
    ttol = 0.1*maxl # 10% is important for curved elements

    if any(X .< Cmin-ttol) || any(X .> Cmax+ttol)
        return false
    end

    # Testing with inverse mapping
    R = inverse_map(shape, C, X, Tol)
    if bdistance(shape, R) > -Tol
        return true
    else
        return false
    end
end

function extrapolator(shape::ShapeType, nips::Int)
    #  Returns a numpy matrix E that extrapolates ip values to nodal values as:
    #
    #                 NodalValues = E * IpValues;
    # where:
    #                             +                +              +
    #                        E = N * (I - EPS1*EPS1 ) + EPS * EPS1
    #
    # and            N = [shape functions matrix]
    #                         1        2        ...        nNodes
    #                  1    [[N_11 N_12
    #                  2     [N_21
    #                  :     [
    #                 nIP    [N_ ...                    ]]
    

    nnodes  = get_nnodes(shape)
    IP      = get_ip_coords(shape, nips)
    ndim    = get_ndim(shape) # shape ndim: not related with the analysis ndim

    #filling N matrix with shape functions of all ips
    N = Array(Float64, nips, nnodes)
    for i=1:nips
        N[i,:] = shape_func(shape, vec(IP[i,:]))
    end

    #calculate extrapolator matrix
    if nips==nnodes
        return inv(N)
    elseif nips>=nnodes
        return pinv(N)
    elseif nips==1
        return pinv(N) # Correction procedure is not applicable for nips==1
    end

    I = eye(nips)

    # εip matrix: Local ip coordinates of integration points
    εip = [ IP[:,1:ndim] zeros(nips) ]

    # ε matrix: Local coordinates of nodal points
    ε = get_local_coords(shape)

    E = pinv(N)*(I - εip*pinv(εip)) + ε*pinv(εip)

    return E
end
