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

Pkg.installed("JSON")==nothing?  Pkg.add("JSON") : nothing
using JSON

export Vect, Matx, vect, matx
export norm2
export subs_equal_by_approx, print_matrix
export RED, GREEN, BOLD, DEFAULT
export @out, @check

# Constants
const RED     = "\x1b[31m"
const GREEN   = "\x1b[32m"
const YELLOW  = "\x1b[33m"
const BLUE    = "\x1b[34m"
const MAGENTA = "\x1b[35m"
const CYAN    = "\x1b[36m"
const WHITE   = "\x1b[37m"
const BOLD    = "\x1b[1m"
const DEFAULT = "\x1b[0m"

@windows_only begin
    const RED     = ""
    const GREEN   = "" 
    const BLUE    = ""
    const MAGENTA = ""
    const CYAN    = ""
    const WHITE   = ""
    const BOLD    = "" 
    const DEFAULT = "" 
end


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

function fix_comparison_scalar(expr::Expr)
    mexpr = copy(expr) # expression to be modified
    tol = 1e-6

    fix_comp = function(expr::Expr)
        symb = expr.args[2]
        a = expr.args[1]
        b = expr.args[3]
        if symb == :(==)
            return :(abs($a-$b) < $tol)
        end
        if symb == :(>=)
            return :($a > $b - $tol)
        end
        if symb == :(>)
            return :($a > $b + $tol)
        end
        if symb == :(<=)
            return :($a < $b + $tol)
        end
        if symb == :(<)
            return :($a < $b - $tol)
        end
        return expr
    end

    if mexpr.head == :comparison
        return fix_comp(mexpr)
    else
        for (i,arg) in enumerate(mexpr.args)
            if typeof(arg)!=Expr; continue end
            if arg.head == :comparison
                mexpr.args[i] = fix_comp(arg)
            else
                mexpr.args[i] = fix_comparison_scalar(arg)
            end
        end
    end
    return mexpr
end

function fix_comparison_arrays(expr::Expr)
    mexpr = copy(expr) # expression to be modified
    tol = 1e-6

    fix_comp = function(expr::Expr)
        symb = expr.args[2]
        a = expr.args[1]
        b = expr.args[3]
        if symb == :(==)
            return :(maximum(abs($a-$b)) < $tol)
        end
        if symb == :(>=)
            return :(minimum($a) > maximum($b) - $tol)
        end
        if symb == :(>)
            return :(minimum($a) > maximum($b) + $tol)
        end
        if symb == :(<=)
            return :(maximum($a) < minimum($b) + $tol)
        end
        if symb == :(<)
            return :(maximum($a) < minimum($b) - $tol)
        end
        return expr
    end

    if mexpr.head == :comparison
        return fix_comp(mexpr)
    else
        for (i,arg) in enumerate(mexpr.args)
            if typeof(arg)!=Expr; continue end
            if arg.head == :comparison
                mexpr.args[i] = fix_comp(arg)
            else
                mexpr.args[i] = fix_comparison_arrays(arg)
            end
        end
    end
    return mexpr
end
