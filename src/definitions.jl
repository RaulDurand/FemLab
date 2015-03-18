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

module Definitions

using JSON

import Base.getindex!

export Vect, Matx, vect, matx
export norm2
export DTable, DBook, push!, getindex, save, loadtable
export subs_equal_by_approx, print_matrix

export @out, @check

# Types
typealias Vect Array{Float64, 1}
typealias Matx Array{Float64, 2}

vect(n::Int64) = Array(Float64,n)
matx(n::Int64, m::Int64) = Array(Float64,n,m)


function find_first_var(exp::Expr)
    for s in exp.args
        if isa(s, Symbol)
            return s
        elseif is(s, Expr)
            return find_first_var(s)
        end
    end
end


# Macro for checkin parameters
macro check(exp, msgs...)
    #msg  = isempty(msgs)? "check failed $exp" : msgs[1]
    var  = find_first_var(exp)
    msg  = isempty(msgs)? "Check failed in expression: $exp" : msgs[1]
    name = string(exp)
    return quote
        if $(esc(exp))
            nothing
        else
            error($(esc(msg)))
        end
    end
end

# Fancy matrix printing
function print_matrix(M::Array{Float64,2})
    n, m = size(M)
    for i=1:n
        for j=1:m
            @printf( "%23.11e", M[i,j] )
        end
        println()
    end
end



# Pseudo determinant of non-square matrices
function norm2(J)

    if ndims(J)==1; return norm(J) end

    r, c = size(J)
    if r==1; return norm(J) end
    if r==2 && c==3
        j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        j2 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
        j3 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
        return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian determinant
    end
    if r==c; return det(J) end
    error("No rule to calculate norm2 of a $r x $c matrix")
end


# DTable object

type DTable
    header::Array{Symbol,1}
    data  ::Array{Array{Float64,1},1}
    dict  ::Dict{Symbol,Array{Float64,1}} # Data index
    function DTable()
        this = new()
        this.header = Array(Symbol,0)
        #this.dict   = []
        return this
    end
    function DTable(header::Array{Symbol}, matrix::Array{Float64}=Float64[])
        this = new()
        this.header = copy(vec(header))
        if length(matrix)==0
            this.data = [ [] for s in header]
        else
            nh = length(header)
            nf = size(matrix,2)
            if nh != nf; error("DTable: header and data fields do not match") end
            this.data = [ matrix[1:end,i] for i=1:nh]
        end
        this.dict = [ k=>v for (k,v) in zip(header, this.data) ]
        return this
    end
end

type DBook
    tables::Array{DTable, 1}
    function DBook()
        this = new()
        this.tables = DTable[]
        return this
    end
end

import Base.push!
function push!(table::DTable, row::Array{Float64})
    for (i,val) in enumerate(row)
        push!(table.data[i], val)
    end
end

function push!(book::DBook, table::DTable)
    push!(book.tables, table)
end

function push!(table::DTable, dict::Dict{Symbol,Float64})
    if length(table.header)==0
        table.header = [ k for k in keys(dict)]
        table.data   = [ [v] for (k,v) in dict ]
        table.dict   = [ k=>v for (k,v) in zip(table.header, table.data) ]
    else
        nrows = length(table.data[1])
        for (k,v) in dict
            # Add data
            if k in table.header
                push!(table[k], v)
            else
                # add new header
                push!(table.header, k)
                new_arr = zeros(nrows)
                push!(table.data, new_arr)
                table.dict[k] = new_arr
                push!(new_arr, v)
            end
        end
        # Add zero for missing values if any
        for arr in table.data
            if length(arr)==nrows
                push!(arr, 0.0)
            end
        end
    end
end

function getindex(table::DTable, field::Symbol)
    return table.dict[field]
end

function getindex(book::DBook, index::Int64)
    return book.tables[index]
end

function save(table::DTable, filename::String; verbose=true, format="dat")
    f  = open(filename, "w")
    nc = length(table.header)   # number of fields (columns)
    nr = length(table.data[1])  # number of rows

    basename, ext = splitext(filename)
    format = (ext == "")? "dat" : ext[2:end]

    if format=="dat"
        # print header
        for i=1:nc
            @printf(f, "%18s", table.header[i])
            print(f, i!=nc? "\t" : "\n")
        end

        # print values
        for i=1:nr
            for j=1:nc
                @printf(f, "%18.10e", table.data[j][i])
                print(f, j!=nc? "\t" : "\n")
            end
        end
        close(f)

        if verbose  println("  file $filename written") end
        return
    end

    if format=="json"
        # generate dictionary
        str  = JSON.json(table.dict, 4)
        print(f, str)
        close(f)

        if verbose  println("  file $filename written (DTable)") end
        return
    end
end

function save(book::DBook, filename::String; verbose=true, format="dat")
    f  = open(filename, "w") 

    if format=="json"
        # generate dictionary
        dict_arr = [ table.dict for table in book.tables ]
        str  = JSON.json(dict_arr, 4)
        print(f, str)
        close(f)

        if verbose  println("  file $filename written (DBook)") end
        return
    end

    if format=="dat" # saves only the last table
        save(book.tables[end], filename, verbose=false)
        if verbose  println("  file $filename written (DBook)") end
        return
    end

end


function loadtable(filename::String)
    data, headstr = readdlm(filename, '\t',header=true)
    header = Symbol[ symbol(strip(field)) for field in headstr ]
    @show typeof(data)

    table = DTable(header, data)
    return table
end


function subs_equal_by_approx(expr::Expr)
    mexpr = copy(expr) # expression to be modified
    for (i,arg) in enumerate(mexpr.args)
        if typeof(arg)!=Expr; continue end
        if arg.head == :comparison
            if arg.args[2] == :(==)
                a = arg.args[1]
                b = arg.args[3]
                mexpr.args[i] = :(isapprox($a,$b,rtol= 1e-8 ))
                continue
            end
        else
            subs_equal_by_approx(arg)
        end
    end
    return mexpr
end


end#module
