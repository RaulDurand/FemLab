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


abstract AbsJoint<:Mechanical

# Return the class of element where this material can be used
client_shape_class(mat::AbsJoint) = JOINT_SHAPE

function init_elem(elem::Element, mat::AbsJoint)

    # Get linked elements
    e1 = elem.linked_elems[1]
    e2 = elem.linked_elems[2]

    # Volume from first linked element
    V1 = 0.0
    C1 = getcoords(e1)
    for ip in e1.ips
        dNdR = e1.shape.deriv(ip.R)
        J    = dNdR*C1
        detJ = det(J)
        V1  += detJ*ip.w
    end

    # Volume from second linked element
    V2 = 0.0
    C2 = getcoords(e2)
    for ip in e2.ips
        dNdR = e2.shape.deriv(ip.R)
        J    = dNdR*C2
        detJ = det(J)
        V2  += detJ*ip.w
    end

    # Area of joint element
    A = 0.0
    C = getcoords(elem)
    n = div(length(elem.nodes), 2)
    C = C[1:n, :]
    bshape = elem.shape.basic_shape
    for ip in elem.ips
        # compute shape Jacobian
        dNdR = bshape.deriv(ip.R)
        J    = dNdR*C
        detJ = norm2(J)

        # compute K
        A += detJ*ip.w
    end

    # Calculate and save h at joint element's integration points
    h = (V1+V2)/(2.0*A)
    for ip in elem.ips
        ip.data.h = h
    end

end

function matrixT(J::Matrix{Float64})
    if size(J,1)==2
        L2 = vec(J[1,:])
        L3 = vec(J[2,:])
        L1 = cross(L2, L3)  # L1 is normal to the first element face
        L1 /= norm(L1)
        L2 /= norm(L2)
        L3 /= norm(L3)
        return [L1 L2 L3]'
    else
        L2 = vec(J)
        L1 = [ L2[2], -L2[1] ] # It follows the anti-clockwise numbering of 2D elements: L1 should be normal to the first element face
        L1 /= norm(L1)
        L2 /= norm(L2)
        return [L1 L2]'
    end
end

function elem_jacobian(mat::AbsJoint, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    bshape = elem.shape.basic_shape

    C = getcoords(elem)[1:hnodes,:]
    B = zeros(ndim, nnodes*ndim)
    K = zeros(nnodes*ndim, nnodes*ndim)

    DB = zeros(ndim, nnodes*ndim)
    J  = zeros(ndim-1, ndim)
    NN = zeros(ndim, nnodes*ndim)

    for ip in elem.ips
        # compute shape Jacobian
        N    = bshape.func(ip.R)
        dNdR = bshape.deriv(ip.R)
        
        @gemm J = dNdR*C
        detJ = norm2(J)

        # compute B matrix
        T    = matrixT(J)
        NN[:,:] = 0.0  # NN = [ -N[]  N[] ]
        for i=1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @gemm B = T*NN

        # compute K
        coef = detJ*ip.w
        D    = mountD(elem.mat, ip.data)
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end

    return K
end

function update!(mat::AbsJoint, elem::Element, dU::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    bshape = elem.shape.basic_shape
    mat    = elem.mat

    dF = zeros(nnodes*ndim)
    C = getcoords(elem)[1:hnodes,:]
    B = zeros(ndim, nnodes*ndim)

    DB = zeros(ndim, nnodes*ndim)
    J  = zeros(ndim-1, ndim)
    NN = zeros(ndim, nnodes*ndim)
    Δω = zeros(ndim)

    for ip in elem.ips
        # compute shape Jacobian
        N    = bshape.func(ip.R)
        dNdR = bshape.deriv(ip.R)
        @gemm J = dNdR*C
        detJ = norm2(J)

        # compute B matrix
        T    = matrixT(J)
        NN[:,:] = 0.0  # NN = [ -N[]  N[] ]
        for i=1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end
        @gemm B = T*NN

        # internal force
        @gemv Δω = B*dU
        Δσ   = stress_update(mat, ip.data, Δω)
        coef = detJ*ip.w
        @gemv dF += coef*B'*Δσ
    end

    return dF
end

function node_and_elem_vals(mat::AbsJoint, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    # Elem vals
    nips = length(elem.ips) 
    E = extrapolator(elem.shape.basic_shape, nips)

    Sn = E*[ ip.data.σ[1] for ip in elem.ips ]
    Wn = E*[ ip.data.w[1] for ip in elem.ips ]

    #Up = [ ip.data.upa for ip in elem.ips ]
    #T  = [ ip.data.σ[2] for ip in elem.ips ]

    node_vals[:sn] = [ Sn; Sn ]
    node_vals[:wn] = [ Wn; Wn ]

    #node_vals[:up] = [ Up; Up ]
    #node_vals[:tau] = [ T; T ]

    #all_ip_vals = [ getvals(mat, ip.data) for ip in elem.ips ]
    #labels      = keys(all_ip_vals[1])
    #nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    #IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    #E = extrapolator(get_basic_shape(shape), nips)
    #N = E*IP # (nnodes x nvals)

    # Filling nodal and elem vals
    #for (i,key) in enumerate(labels)
        #node_vals[key] = N[:,i]
    #end

    return node_vals, elem_vals

end

