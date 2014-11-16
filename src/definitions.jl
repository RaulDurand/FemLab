module Definitions

export Vect, Matx, vect, matx
export norm2
export DTable, push!, getindex, save, loadtable

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
    msg  = isempty(msgs)? "Invalid value for $var since it does not satisfy $exp" : msgs[1]
    name = string(exp)
    return quote
        if $(esc(exp))
            nothing
        else
            error($(esc(msg)))
        end
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
    dict  ::Dict{Symbol,Int} # Data index
    function DTable()
        this = new()
        this.header = Array(Symbol,0)
        return this
    end
    function DTable(header::Array{Symbol}, data::Array{Float64}=Float64[])
        this = new()
        this.header = copy(vec(header))
        if length(data)==0
            this.data = [ [] for s in header]
        else
            nh = length(header)
            nf = size(data,2)
            if nh != nf; error("DTable: header and data fields do not match") end
            this.data = [ data[1:end,i] for i=1:nh]
        end
        this.dict   = [ s=>i for (i,s) in enumerate(header) ]
        return this
    end
end

import Base.push!
function push!(table::DTable, row::Array{Float64})
    for (i,val) in enumerate(row)
        push!(table.data[i], val)
    end
end

function push!(table::DTable, dict::Dict{Symbol,Float64})
    if length(table.header)==0
        table.header = [ k for k in keys(dict)]
        table.data   = [ [v] for (k,v) in dict ]
        table.dict    = [ s=>i for (i,s) in enumerate(table.header) ]
    else
        for (k,v) in dict
            push!(table[k], v)
        end
    end
end

import Base.getindex!
function getindex(table::DTable, field::Symbol)
    index = table.dict[field]
    return table.data[index]
end

function save(table::DTable, filename::String, verbose=false)
    f = open(filename, "w")
    nc = length(table.header)
    nr = length(table.data[1])

    # print header
    for i=1:nc
        @printf(f, "%-18s", table.header[i])
        print(f, i!=nc? "\t" : "\n")
    end

    # print values
    nc = length(table.data)
    #nr, nc = size(table.data)
    for i=1:nr
        for j=1:nc
            @printf(f, "%18.10e", table.data[j][i])
            print(f, j!=nc? "\t" : "\n")
        end
    end
    close(f)

    if verbose
        println("  file $filename written")
    end
end


function loadtable(filename::String)
    data, headstr = readdlm(filename, '\t',header=true)
    header = Symbol[ symbol(strip(field)) for field in headstr ]

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
                mexpr.args[i] = :(isapprox($a,$b,rtol=1.e-8))
                continue
            end
        else
            subs_equal_by_approx(arg)
        end
    end
    return mexpr
end


end#module
