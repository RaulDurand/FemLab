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

abstract BC

type NodeBC <: BC
    expr ::Expr
    conds::Array
    nodes::Array{Node,1}

    function NodeBC(expr::Expr; conds::Array...)
        this = new()
        this.expr = expr
        this.conds = conds
        this.nodes = []
        return this
    end

    function NodeBC(nodes::Array{Node,1}; conds::Any...)
        this = new()
        this.expr = :()
        this.nodes = nodes
        this.conds = conds
        return this
    end

    function NodeBC(node::Node; conds::Any...)
        return NodeBC( [node]; conds... )
    end

end

type FaceBC <: BC
    expr ::Expr
    conds::Array
    faces::Array{Face,1}

    function FaceBC(expr::Expr; conds::Array...)
        this = new()
        this.expr = expr
        this.conds = conds
        return this
    end

    function FaceBC(faces::Array{Face,1}; conds::Any...)
        this = new()
        this.expr = :()
        this.faces = faces
        this.conds = conds
        return this
    end

    function FaceBC(face::Face; conds::Any...)
        return FaceBC( [face]; conds... )
    end
end

type EdgeBC <: BC
    expr ::Expr
    conds::Array
    edges::Array{Edge,1}

    function EdgeBC(expr::Expr; conds::Array...)
        this = new()
        this.expr = expr
        this.conds = conds
        return this
    end

    function EdgeBC(edges::Array{Edge,1}; conds::Any...)
        this = new()
        this.expr = :()
        this.edges = edges
        this.conds = conds
        return this
    end

    function EdgeBC(edge::Edge; conds::Any...)
        return EdgeBC( [edge]; conds... )
    end
end



# Define boundary conditions for a node
"""
`apply_bc(node, bcs...)`

Applies one or several boundary conditions `bcs` to a `node` object.
In a mechanical analysis, essential and natural boundary conditions can be set using this function.
For example:
```
node = Node([1.0, 1.0, 1.0])
apply_bc(node, fx=10.0, uy=0.01)

```
"""
function apply_bc(node::Node; args...)
    for (key,val) in args
        if !haskey(node.dofdict, key); error("key ($key) not found in node ($(node.id)).") end
        dof = node.dofdict[key]
        if key==dof.sU
            dof.prescU = true
            dof.bryU  = val
        else
            dof.bryF += val
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

# Define boundary conditions for a collection of nodes
"""
`apply_bc(nodes, bcs...)`

Applies one or several boundary conditions `bcs` to a set of Node objects `nodes`.
"""
function apply_bc(nodes::Array{Node,1}; args...)
    length(nodes)==0 && warn("Warning, applying boundary conditions to empty array of nodes\n")

    for node in nodes
        apply_bc(node; args...)
    end
end

# Define boundary conditions for a face
function apply_bc(face::Face; args...)
    oelem = face.oelem # owner element
    oelem==nothing && error("Face with no owner element")

    for (key,val) in args
        apply_facet_bc(oelem.mat, oelem, face, key, float(val))
    end
end

# Define boundary conditions for a collection of faces
function apply_bc(faces::Array{Face,1}; args...)
    length(faces)==0 && warn("Warning, applying boundary conditions to empty array of faces\n")

    for face in faces
        apply_bc(face; args...)
    end
end
