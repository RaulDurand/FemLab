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

export set_mat

abstract Mechanical<:Material

function config_dofs(::Mechanical, elem::Element)
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        if elem.ndim==3; add_dof(node, :uz, :fz) end
    end
end

function get_map(elem::Element)
    get_map(elem.mat, elem)
end

function get_map(::Mechanical, elem::Element)
    dof_keys = (:ux, :uy, :uz)[1:elem.ndim]
    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)
end

function set_facet_bc(mat::Mechanical, oelem::Element, face::Face, key::Symbol, val::Float64)
    fnodes = face.nodes
    fshape = face.shape
    ndim   = oelem.ndim

    if key in (:fx, :fy, :fz)
        error("Boundary condition $key is not allowed in a facet; consider tx, ty, tz or tn")
    end

    if !(key in (:ux, :uy, :uz, :tx, :ty, :tz, :tn))
        error("Boundary condition $key is not applicable in a facet of an element with material $(typeof(mat))")
    end

    if (key in (:tz, :uz) && ndim==2)
        error("Boundary condition $key is not applicable in a 2D analysis")
    end

    # Apply the boundary conditions
    if key in (:ux, :uy, :uz)
        for node in fnodes
            dof = node.dofdict[key]
            dof.prescU = true
            dof.bryU   = val
        end
        return
    end

    # Force boundary condition
    nfnodes = length(fnodes)

    # Calculate the face coordinates matrix
    C = getcoords(fnodes, ndim)

    # Calculate the vector with values to apply
    V = zeros(ndim)
    if key == :tx ; V[1] = val end
    if key == :ty ; V[2] = val end
    if key == :tz ; V[3] = val end

    # Calculate the nodal values
    F = zeros(nfnodes, ndim)

    f_ips = get_ip_coords(fshape)

    for i=1:size(f_ips,1)
        R = vec(f_ips[i,:])
        w = R[end]
        N = shape_func(fshape, R)
        D = deriv_func(fshape, R)
        J = D*C
        nJ = norm2(J)

        if key == :tn && ndim==2
            n = [J[1,2], -J[1,1]]
            V = val*n/norm(n)
        end
        if key == :tn && ndim==3
            n = cross(J[1,:], J[2,:])
            V = val*n/norm(n)
        end

        F += N*V'*(nJ*w)
    end

    # Setting bc into nodes
    for (i,node) in enumerate(fnodes)
        node.dofdict[:fx].bryF += F[i,1]
        node.dofdict[:fy].bryF += F[i,2]
        if ndim==3; node.dofdict[:fz].bryF += F[i,3] end
    end
end

function init_elem(elem::Element, mat::Mechanical)
    # empty, called by set_mat(...)
end

function setB(ndim::Int, dNdX::Matx, detJ::Float64, B::Matx)
    nnodes = size(dNdX,2)
    sqr2 = √2.0
    B[:] = 0.0

    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[1,i]
            B[2,2+j*ndim] = dNdX[2,i]
            B[4,1+j*ndim] = dNdX[2,i]/sqr2; B[4,2+j*ndim] = dNdX[1,i]/sqr2
        end
        if gl_stress_state==:axisymmetric # TODO: Check this
            for i in 1:nnodes
                N = shape_func(elem.shape, R)
                j = i-1
                r = R[0]
                B[1,1+j*ndim] = dNdX[1,i]
                B[2,2+j*ndim] = dNdX[2,i]
                B[3,1+j*ndim] =    N[i]/r
                B[4,1+j*ndim] = dNdX[2,i]/sqr2; B[4,2+j*ndim] = dNdX[1,i]/sqr2
            end
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[1,i]
            dNdy = dNdX[2,i]
            dNdz = dNdX[3,i]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,1+j*ndim] = dNdy/sqr2;   B[4,2+j*ndim] = dNdx/sqr2
            B[5,2+j*ndim] = dNdz/sqr2;   B[5,3+j*ndim] = dNdy/sqr2
            B[6,1+j*ndim] = dNdz/sqr2;   B[6,3+j*ndim] = dNdx/sqr2
        end
    end

    return detJ
end

function elem_jacobian(::Mechanical, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    mat    = elem.mat
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)

    DB = Array(Float64, 6, nnodes*ndim)
    J  = Array(Float64, ndim, ndim)
    dNdX = Array(Float64, ndim, nnodes)

    for ip in elem.ips

        # compute B matrix
        dNdR = deriv_func(elem.shape, ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        setB(ndim, dNdX, detJ, B)

        # compute K
        coef = detJ*ip.w
        D    = calcD(mat, ip.data) 
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end
    return K
end

function elem_jacobian(elem::Element)
    elem_jacobian(elem.mat, elem)
end

function elem_RHS(elem::Element)
    return elem_RHS(elem.mat, elem)
end

function elem_RHS(::Mechanical, elem::Element)
    # Appropriate for body forces
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    return zeros(ndim*nnodes)
end

function update!(elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    update!(elem.mat, elem, DU, DF)
end

function update!(::Mechanical, elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    map    = get_map(elem)
    mat    = elem.mat

    ΔF = zeros(nnodes*ndim)
    ΔU = DU[map]
    B  = zeros(6, nnodes*ndim)

    DB = Array(Float64, 6, nnodes*ndim)
    J  = Array(Float64, ndim, ndim)
    dNdX = Array(Float64, ndim, nnodes)
    Δε = zeros(6)

    C = getcoords(elem)
    for ip in elem.ips

        # compute B matrix
        dNdR = deriv_func(elem.shape, ip.R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        setB(ndim, dNdX, detJ, B)

        @gemv Δε = B*ΔU
        Δσ   = stress_update(mat, ip.data, Δε)
        coef = detJ*ip.w
        @gemv ΔF += coef*B'*Δσ
    end

    # Update global vector
    DF[map] += ΔF

end



function elem_vals(mat::Mechanical, elem::Element)
    return Dict{Symbol, Float64}()
end

function node_and_elem_vals(mat::Mechanical, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    e_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    # values from integration points
    all_ip_vals = [ getvals(mat, ip.data) for ip in elem.ips ]
    labels      = keys(all_ip_vals[1])
    nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    # extrapolate values to nodes
    E = extrapolator(elem.shape, nips)
    N = E*IP # (nnodes x nvals)


    # Filling nodal and elem vals
    for (i,key) in enumerate(labels)
        node_vals[key] = N[:,i]
    end

    e_vals = elem_vals(mat, elem)

    return node_vals, e_vals

end

include("tensors.jl")

include("solid.jl")
include("dp.jl")
include("kotsovos.jl")
include("mazars.jl")

include("truss.jl")
include("pptruss.jl")

include("joint.jl")
include("hordijk.jl")
include("nljoint.jl")

include("joint1d.jl")
include("mcjoint1d.jl")
include("cebjoint1d.jl")

#include("bilinear.jl")
#include("smeared_crack.jl")
include("mcjoint.jl")
