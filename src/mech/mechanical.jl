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

abstract Mechanical<:Material

function config_dofs(::Mechanical, elem::Element)
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        if elem.ndim==3; add_dof(node, :uz, :fz) end
    end
end

function init_elem(elem::Element, mat::Mechanical)
    # empty, called by set_mat(...)
    # to be speciallized according to materials
end

function get_map(elem::Element)
    get_map(elem.mat, elem)
end

function get_map(::Mechanical, elem::Element)
    dof_keys = (:ux, :uy, :uz)[1:elem.ndim]
    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)
end

function elem_jacobian(elem::Element)
    elem_jacobian(elem.mat, elem)
end

function elem_RHS(elem::Element)
    return elem_RHS(elem.mat, elem)
end

function elem_RHS(::Mechanical, elem::Element)
    # Appropriate for body forces
    return zeros(length(get_map(elem)))
end

function update!(elem::Element, dU::Array{Float64,1})
    return update!(elem.mat, elem, dU)
end

function elem_vals(mat::Mechanical, elem::Element)
    return Dict{Symbol, Float64}()
end

