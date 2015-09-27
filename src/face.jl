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

export Face
export @get_faces


type Face
    shape::ShapeType
    nodes::Array{Node,1}
    ndim ::Integer
    oelem::Union(Element,Nothing)
    isedge::Bool
    function Face(shape, nodes, ndim)
        this = new(shape, nodes, ndim)
        this.oelem  = nothing
        this.isedge = false
        return this
    end
end

typealias Edge Face

function getindex(faces::Array{Face,1}, cond::Expr)
    condm = fix_comparison_arrays(cond)
   
    funex = :( (x,y,z) -> x*y*z )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun   = eval(funex)
    catch
        error("Faces getindex: Invalid condition ", cond)
    end

    result = Array(Face,0)
    for face in faces
        coords = getcoords(face.nodes)
        x = coords[:,1]
        y = coords[:,2]
        z = coords[:,3]
        if fun(x, y, z)
            push!(result, face) 
        end
    end
    return result
end

getindex(faces::Array{Face,1}, cond::String) = getindex(faces, parse(cond))

function set_bc(face::Face; args...)
    oelem = face.oelem # owner element
    if oelem==nothing; error("Face with no owner element") end

    for (key,val) in args
        set_facet_bc(oelem.mat, oelem, face, key, float(val))
    end
end

function set_bc(faces::Array{Face,1}; args...)
    if length(faces)==0; pcolor(:red, "Warning, applying boundary conditions to empty array of faces\n") end
    for face in faces
        set_bc(face; args...)
    end
end

# Macro to filter faces using a condition expression
macro get_faces(dom, expr)

    # fix condition
    cond = fix_comparison_arrays(expr)

    # generate the filter function
    func = quote
        (f::Face) -> begin 
            x = [ node.X[1] for node in f.nodes ]
            y = [ node.X[2] for node in f.nodes ]
            z = [ node.X[3] for node in f.nodes ]
            $(cond) 
        end
    end

    quote
        ff = $(esc(func))
        tt = Bool[ ff(f) for f in $(esc(dom)).faces]
        $(esc(dom)).faces[ tt ]
    end
end
