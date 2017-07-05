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


abstract type AbsJoint1D<:Mechanical end

# Return the class of element where this material can be used
client_shape_class(mat::AbsJoint1D) = JOINT1D_SHAPE


function mount_T(J::Matx)
    ndim = length(J)
    nJ   = norm(J)
    L1   = vec(J/nJ)

    if ndim==2
        L2 = [ -L1[2],  L1[1] ]
        return hcat(L1,L2)' # TODO: check for 2D
    end

    # Finding second vector
    if     abs(L1[1]) == 1.0; L2 = [0.0, 1.0, 0.0]
    elseif abs(L1[2]) == 1.0; L2 = [0.0, 0.0, 1.0]
    elseif abs(L1[3]) == 1.0; L2 = [1.0, 0.0, 0.0]
    else
        # Auxiliar vector L which must be different from L1
        L = [1.0, 0.0, 0.0]
        if norm(L-L1) < 1.0e-4; L = [0.0, 1.0, 0.0] end
        # Performing cross product to obtain a second vector
        L2  = cross(L1, L)
        L2 /= norm(L2)
    end

    # Finding third vector
    L3 = cross(L1, L2)
    L3 /= norm(L3)

    return hcat(L1, L2, L3)'
end


function mountB(mat::AbsJoint1D, elem::Element, R, Ch, Ct)
    # Calculates the matrix that relates nodal displacements with relative displacements
    # ==================================================================================

    # B = T* [NN*MM  -NN]      ndim x ndim*(m+n)

    # where
    # T is a direction cosines matrix
    # NN is a matrix containing truss element shape functions
    # evaluated at the point of interest R.
    # MM is a matrix containing tresspased element shape functions
    # evaluated at the n truss nodal points.

    #          [ M_11*I M_21*I ... M_m1*I]
    #     MM = [ M_12*I M_22*I ... M_m2*I]
    #          [ M_13*I M_23*I ... M_mn*I] ndim*n x ndim*m

    # where
    # M_12 is the first shape function from the tresspased solid element
    # evaluated at the second node of the truss element.
    # I is a ndim x ndim identity matrix


    ndim = elem.ndim
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    nnodes  = length(elem.nodes)
    nbnodes = length(bar.nodes)
    nsnodes = length(hook.nodes)
    D = bar.shape.deriv(R)
    J = D*Ct
    T = mount_T(J)

    # Mount NN matrix
    N = bar.shape.func(R)
    NN = hcat([ Ni*eye(ndim) for Ni in N  ]...)

    # Mount MM matrix
    stack = Array{Float64,2}[]
    for i=1:nbnodes
        Xj = bar.nodes[i].X
        R  = inverse_map(hook.shape, Ch, Xj)
        M  = hook.shape.func(R)
        for Mi in M
            push!(stack, Mi*eye(ndim))
        end
    end
    MM = hvcat(nsnodes, stack...)

    B = -T*[ NN*MM  -NN ]
    detJ = norm(J)
    return B, detJ
end

function elem_init(mat::AbsJoint1D, elem::Element)::Void
    B_detJ = []
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    Ch = elem_coords(hook)
    Ct = elem_coords(bar)
    for ip in elem.ips
        B, detJ = mountB(elem.mat, elem, ip.R, Ch, Ct)
        push!(B_detJ, (B, detJ))
    end
    elem.cache[:B_detJ] = B_detJ
    return nothing
end

function elem_stiffness(mat::AbsJoint1D, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    #Ch = elem_coords(hook)
    #Ct = elem_coords(bar)

    K  = zeros(nnodes*ndim, nnodes*ndim)
    DB = zeros(ndim, nnodes*ndim)

    for (i,ip) in enumerate(elem.ips)
        #B, detJ = mountB(elem.mat, elem, ip.R, Ch, Ct)
        B, detJ = elem.cache[:B_detJ][i]
        D    = calcD(elem.mat, ip.data)
        D[1,1]*=mat.h
        coef = detJ*ip.w
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end
    return K
end

function elem_dF!(mat::AbsJoint1D, elem::Element, dU::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    mat    = elem.mat

    dF = zeros(nnodes*ndim)
    #B  = zeros(ndim, nnodes*ndim)
    deps = zeros(ndim)

    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    #Ct   = elem_coords(bar)
    #Ch   = elem_coords(hook)
    #for ip in elem.ips
    for (i,ip) in enumerate(elem.ips)
        #detJ = mountB(elem.mat, elem, ip.R, Ch, Ct, B)
        B, detJ = elem.cache[:B_detJ][i]
        D    = calcD(mat, ip.data)
        @gemv deps = B*dU
        dsig = stress_update(mat, ip.data, deps)
        coef = detJ*ip.w
        dsig[1]  *= mat.h
        @gemv dF += coef*B'*dsig
    end

    return dF
end

function elem_and_node_vals(mat::AbsJoint1D, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    # Elem vals
    all_ip_vals = [ ip_state_vals(mat, ip.data) for ip in elem.ips ]
    labels      = keys(all_ip_vals[1])
    nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]

    E = extrapolator(bar.shape, nips)
    N = E*IP # (nnodes x nvals)

    nhnodes = length(hook.nodes)
    N = vcat( zeros(nhnodes, length(labels)), N )

    # Filling nodal and elem vals
    for (i,key) in enumerate(labels)
        node_vals[key] = N[:,i]
    end

    return node_vals, elem_vals

end

