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


export EmbTruss

abstract AbsEmbTruss<:Mechanical

# Return the class of element where this material can be used
client_elem_class(mat::AbsEmbTruss) = :EMBEDDED

function config_dofs(::AbsEmbTruss, elem::Element)
    # No-op function.
    # The solid linked element will set the required dofs.
end

function get_map(::AbsEmbTruss, elem::Element)
    # Return the map from the solid linked element
    return get_map(elem.linked_elems[1])
end

function mountNN(mat::AbsEmbTruss, elem::Element)
    solid = elem.linked_elems[1]
    ndim = elem.ndim
    n  = length(solid.nodes)
    m  = length(elem.nodes)
    NN = zeros(ndim*n, ndim*m)
    Cs = getcoords(solid)

    for j=1:m 
        R = inverse_map(solid.shape, Cs, elem.nodes[j].X)
        N = shape_func(solid.shape, R)
        for i=1:n
            for k=1:ndim
                NN[(i-1)*ndim+k, (j-1)*ndim+k] = N[i]
            end
        end
    end
    return NN
end

function elem_jacobian(mat::AbsEmbTruss, elem::Element)
    Kb = elem_jacobian(mat.trussmat, elem)
    NN = mountNN(mat, elem)
    return NN*Kb*NN' 
end

function update!(mat::AbsEmbTruss, elem::Element, dU::Array{Float64,1})
    NN   = mountNN(mat, elem)
    dUtr = NN'*dU
    dFtr = update!(mat.trussmat, elem, dUtr)
    dF   = NN*dFtr
    return dF
end

function node_and_elem_vals(mat::AbsEmbTruss, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    solid = elem.linked_elems[1]

    # get solid nodal displacements
    n = length(solid.nodes)
    dU = zeros(ndim*n)
    for (i,node) in enumerate(solid.nodes)
        for (j,dof) in enumerate(node.dofs)
            dU[(i-1)*ndim+j] = dof.U
        end
    end

    # get nodal displacements at embedded truss nodes by interpolation
    NN = mountNN(mat, elem)
    dUtr = NN'*dU

    # Fill nodal values
    for (i,key) in enumerate( (:ux, :uy, :uz)[1:ndim] )
        node_vals[key] = dUtr[i:ndim:end]
    end

    # Elem vals
    all_ip_vals = [ getvals(mat, ip.data) for ip in elem.ips ]

    # completing with axial forces
    A = mat.trussmat.A
    for ip_val in all_ip_vals
        ip_val[:A ] = A
        ip_val[:Fa] = A*ip_val[:sa]
    end
    labels      = keys(all_ip_vals[1])
    nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    E = extrapolator(elem.shape, nips)
    N = E*IP # (nnodes x nvals)

    # Filling nodal and elem vals
    for (i,key) in enumerate(labels)
        #node_vals[key] = N[:,i] # No nodal values for truss
        elem_vals[key] = mean(IP[:,i])
    end

    return node_vals, elem_vals
end
