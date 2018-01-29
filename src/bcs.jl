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

abstract type BC end


###############################################################################
# Concrete structures for Boudary Conditions


mutable struct NodesBC <: BC
    expr ::Expr
    conds::Array
    nodes::Array{Node,1}

    function NodesBC(expr::Expr; conds::Array...)
        return new(expr, conds, [])
    end

    function NodesBC(nodes::Array{Node,1}; conds::Any...)
        return new(:(), conds, nodes)
    end

    function NodesBC(node::Node; conds::Any...)
        return NodesBC( [node]; conds... )
    end

end


mutable struct ElementBC <: BC
  
    expr ::Expr
    conds::Array
    elems::Array{Element,1}

    function ElementBC(expr::Expr; conds::Array...)
        return new(expr, conds, [])
    end

    function ElementBC(elems::Array{Element,1}; conds::Any...)
        return new(:(), conds, elems)
    end

    function ElementBC(elem::Element; conds::Any...)
        return NodesBC( [elem]; conds... )
    end

end


abstract type FacetsBC<:BC end

mutable struct FacesBC<:FacetsBC
    expr  ::Expr
    conds ::Array
    faces::Array{Face,1}

    function FacesBC(expr::Expr; conds::Array...)
        return new(expr, conds, [])
    end

    function FacesBC(faces::Array{Face,1}; conds::Any...)
        return new(:(), conds, faces)
    end

    function FacesBC(face::Face; conds::Any...)
        return new(:(), conds, [face])
    end
end

mutable struct EdgesBC<:FacetsBC
    expr  ::Expr
    conds ::Array
    edges::Array{Edge,1}

    function EdgesBC(expr::Expr; conds::Array...)
        return new(expr, conds, [])
    end

    function EdgesBC(edges::Array{Edge,1}; conds::Any...)
        return new(:(), conds, edges)
    end

    function EdgesBC(edge::Edge; conds::Any...)
        return new(:(), conds, [edge])
    end
end


# For compatibility with older versions
NodeBC = NodesBC
FaceBC = FacesBC
EdgeBC = EdgesBC


###############################################################################
# Functions to setup boundary conditions

function setup_bc!(bc::NodesBC, domain)
    if bc.expr != :(); bc.nodes = domain.nodes[bc.expr] end
end

function setup_bc!(bc::FacesBC, domain)
    if bc.expr != :(); bc.faces = domain.faces[bc.expr] end
end

function setup_bc!(bc::EdgesBC, domain)
    if bc.expr != :(); bc.edges = domain.edges[bc.expr] end
end


###############################################################################
# Functions to apply boundary conditions


"""
`apply_bc(bc::NodesBC)`

Applies one or several boundary conditions `bcs` to a `node` object.
In a mechanical analysis, essential and natural boundary conditions can be set using this function.
For example:
```
node = Node([1.0, 1.0, 1.0])
apply_bc(node, fx=10.0, uy=0.01)

```
"""
function apply_bc(bc::NodesBC)
    length(bc.nodes)==0 && warn("Warning, applying boundary conditions to empty array of nodes")
    for node in bc.nodes
        for (key,val) in bc.conds
            if !haskey(node.dofdict, key); error("key ($key) not found in node ($(node.id)).") end
            dof = node.dofdict[key]
            if key==dof.sU
                dof.prescU = true
                dof.bryU  = val
            else
                #dof.bryF += val
                dof.bryF = Expr(:call, :+, dof.bryF, val)
            end
        end
    end
end


function apply_bc{T<:FacetsBC}(bc::T)
    facets = T==FacesBC ? bc.faces : bc.edges
    length(facets)==0 && warn("$T: applying boundary conditions to empty array of facets")
    for facet in facets
        oelem = facet.oelem # owner element
        oelem==nothing && error("Facet with no owner element")

        for (key,val) in bc.conds
            apply_facet_bc(oelem.mat, oelem, facet, key, val)
        end
    end
end

function apply_bc(bc::ElementBC)
    
    for elem in bc.elems
        for (key,val) in bc.conds
            apply_elem_bc(elem.mat, elem, key, float(val))
        end
    end
end


# Clear all boundary conditions in a node
"""
`clear_bc(node)`

Clears all boundary conditions previously set in a Node object.
"""
function clear_bc(node::Node)
    for dof in node.dofs
        dof.bryU   = 0.0
        dof.bryF   = 0.0
        dof.prescU = false
        dof.eq_id  = -1
    end
end

