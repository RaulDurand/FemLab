
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

import Base.getindex
import Base.keys
export DTable, DBook, push!, keys, getindex, save, loadtable


# DTable object

type DTable
    data  ::Array{Array{Float64,1},1}
    dict  ::Dict{Symbol,Array{Float64,1}} # Data index
    function DTable()
        this = new()
        this.dict = Dict{Symbol,Array{Float64,1}}() 
        return this
    end
    function DTable(header::Array{Symbol,1}, matrix::Array{Float64}=Float64[])
        this = new()
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

function keys(table::DTable)
    return keys(table.dict)
end

function push!(table::DTable, dict::Dict{Symbol,Float64})
    if length(table.dict)==0
        table.data   = [ [v] for (k,v) in dict ]
        table.dict   = [ k=>v for (k,v) in zip(keys(dict), table.data) ]
    else
        nrows = length(table.data[1])
        for (k,v) in dict
            # Add data
            if haskey(table.dict, k)
                push!(table[k], v)
            else
                # add new column
                new_arr = zeros(nrows)
                push!(new_arr, v)
                push!(table.data, new_arr)
                table.dict[k] = new_arr
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


function save(table::DTable, filename::AbstractString; verbose=true, format="")
    f  = open(filename, "w")
    nc = length(table.dict)   # number of fields (columns)
    nr = length(table.data[1])  # number of rows

    basename, ext = splitext(filename)
    if format==""
        format = (ext == "")? "dat" : ext[2:end]
    end

    if format=="dat"
        # map for ordered header
        symbols = collect(keys(table.dict))
        dmap = [ tup[2] for tup in sort(collect(zip(symbols, 1:nc))) ]

        # print header
        for i=1:nc
            @printf(f, "%18s", symbols[dmap[i]])
            print(f, i!=nc? "\t" : "\n")
        end

        # print values
        for i=1:nr
            for j=1:nc
                @printf(f, "%18.10e", table.data[dmap[j]][i])
                print(f, j!=nc? "\t" : "\n")
            end
        end
        close(f)

        if verbose  pcolor(:green, "  file $filename written\n") end
        return
    end

    if format=="json"
        # enconding
        str  = JSON.json(table.dict, 4)
        print(f, str)
        close(f)

        if verbose  pcolor(:green, "  file $filename written (DTable)\n") end
        return
    end
end


function save(book::DBook, filename::AbstractString; verbose=true, format="dat")
    f  = open(filename, "w") 

    if format=="json"
        # generate dictionary
        dict_arr = [ table.dict for table in book.tables ]
        str  = JSON.json(dict_arr, 4)
        print(f, str)
        close(f)

        if verbose  pcolor(:green, "  file $filename written (DBook)\n") end
        return
    end

    if format=="dat" # saves only the last table
        #save(book.tables[end], filename, verbose=false)
        #if verbose  pcolor(:green, "  file $filename written (DBook)\n") end

        f  = open(filename, "w")

        basename, ext = splitext(filename)
        if format==""
            format = (ext == "")? "dat" : ext[2:end]
        end

        for table in book.tables
            symbols = collect(keys(table.dict))
            nc = length(table.dict)     # number of fields (columns)
            nr = length(table.data[1])  # number of rows

            # map for ordered header
            dmap = [ tup[2] for tup in sort(collect(zip(symbols, 1:nc))) ]

            # print header
            for i=1:nc
                @printf(f, "%18s", symbols[dmap[i]])
                print(f, i!=nc? "\t" : "\n")
            end

            # print values
            for i=1:nr
                for j=1:nc
                    @printf(f, "%18.10e", table.data[dmap[j]][i])
                    print(f, j!=nc? "\t" : "\n")
                end
            end
            print(f, "\n")
        end

        close(f)
        if verbose  pcolor(:green, "  file $filename written\n") end

        return
    end

end


function loadtable(filename::AbstractString)
    data, headstr = readdlm(filename, '\t',header=true)
    header = Symbol[ symbol(strip(field)) for field in headstr ]

    table = DTable(header, data)
    return table
end
