
function subs_equal_by_approx(expr::Expr) # deprecated
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

# Fixes comparisons expressions using a tolerance
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

# Fixes comparisons in array expressions using a tolerance
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
