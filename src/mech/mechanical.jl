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

abstract type Mechanical<:Material end

"""
`elem_config_dofs(mat, elem)`

Sets up the dofs for all nodes in `elem` according to material mat.
This function can be overloaded by concrete types.
"""
function elem_config_dofs(mat::Mechanical, elem::Element)::Void
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        if elem.ndim==3; add_dof(node, :uz, :fz) end
    end
end

"""
`elem_init(mat, elem)`

Sets up `elem` according to material `mat`.
This function is called after mat is assigned to `elem` by function `set_mat`.
This function can be overloaded by concrete types.
"""
function elem_init(mat::Material, elem::Element)::Void
    # No-op function but can be specialized by concrete types
    # This is called by set_mat(...) function
    return nothing
end

"""
`elem_map(mat, elem)`

Returns an array of corresponding equation numbers for all dofs in `elem` according to `mat`.
This function can be overloaded by concrete types.
"""
function elem_map(mat::Mechanical, elem::Element)::Array{Int,1}
    dof_keys = (:ux, :uy, :uz)[1:elem.ndim]
    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)
end

"""
`elem_stiffness(mat, elem)`

Returns the stiffness matrix for `elem` according to material `mat`.
This function must be redefined by concrete types.
"""
function elem_stiffness(mat::Mechanical, elem::Element)::Array{Float64,2}
    error("elem_stiffness function not defined for material type $(typeof(mat))")
end

"""
`elem_dF!(mat, elem)`

Returns the force increment vector dF given a displecement increment vector `dU`
for `elem` according to material `mat`.
This function also updates strains, stresses and internal variables of all
`IpState` objects at integration points.
This function must be redefined by concrete types.
"""
function elem_dF!(mat::Mechanical, elem::Element, dU::Array{Float64,1})
    error("elem_dF function not defined for material type $(typeof(mat))")
end

"""
`elem_RHS(mat, elem)`

Returns the right-hand-side vector for `elem`.
This function can be overloaded by concrete types.
"""
function elem_RHS(mat::Mechanical, elem::Element)::Array{Float64,1}
    # Appropriate for body forces
    return zeros(length(elem_map(mat, elem)))
end

"""
`elem_vals(mat, elem)`

Returns two dictionaries with keys and values for all nodes and for `elem` itself.
Values for the element are intended to be constant along the element.
This function can be overloaded by concrete types.
"""
function elem_and_node_vals(mat::Mechanical, elem::Element)::( Dict{Symbol, Float64}, Dict{Symbol, Float64})
    return Dict{Symbol, Float64}, Dict{Symbol, Float64}()
end
