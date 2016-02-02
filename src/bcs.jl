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

export NodeBC, FaceBC, EdgeBC

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
