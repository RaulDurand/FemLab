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


abstract type AbsEmbTruss<:Mechanical end

# Return the class of element where this material can be used
client_shape_class(mat::AbsEmbTruss) = EMBEDDED

function elem_config_dofs(::AbsEmbTruss, elem::Element)
    # No-op function.
    # The solid linked element will set the required dofs.
end

function elem_map(mat::AbsEmbTruss, elem::Element)
    # Return the map from the solid linked element
    linked = elem.linked_elems[1]
    return elem_map(linked.mat, linked)  #TODO: solve the case when the linked elem has more dofs e.g. temperature.
end

function mountNN(mat::AbsEmbTruss, elem::Element)
    solid = elem.linked_elems[1]
    ndim = elem.ndim
    n  = length(solid.nodes)
    m  = length(elem.nodes)
    NN = zeros(ndim*n, ndim*m)
    Cs = elem_coords(solid)

    for j=1:m 
        R = inverse_map(solid.shape, Cs, elem.nodes[j].X)
        N = solid.shape.func(R)
        for i=1:n
            for k=1:ndim
                NN[(i-1)*ndim+k, (j-1)*ndim+k] = N[i]
            end
        end
    end
    return NN
end


function elem_init(mat::AbsEmbTruss, elem::Element)::Void
    elem.cache[:NN] = mountNN(mat, elem)
    return nothing
end


function elem_stiffness(mat::AbsEmbTruss, elem::Element)
    Kb = elem_stiffness(mat.trussmat, elem)
    NN = elem.cache[:NN]
    return NN*Kb*NN' 
end

function elem_dF!(mat::AbsEmbTruss, elem::Element, dU::Array{Float64,1})
    NN = elem.cache[:NN]
    dUtr = NN'*dU
    dFtr = elem_dF!(mat.trussmat, elem, dUtr)
    dF   = NN*dFtr
    return dF
end

function elem_and_node_vals(mat::AbsEmbTruss, elem::Element)
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
    all_ip_vals = [ ip_state_vals(mat, ip.data) for ip in elem.ips ]

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
