export AbsBeam
abstract type AbsBeam<:Mechanical end

function beam_shape_func(ξ::Float64, nnodes::Int)
    if nnodes==2
        N = Array{Float64}(4)
        x = (ξ+1)/2
        N[1] = 1 - 3*x^2 + 2*x^3
        N[2] = x - 2*x^2 + x^3
        N[3] = 3*x^2 - 2*x^3
        N[4] = x^3 - x^2
    else
        N = Array{Float64}(4)
    end
    return N
end

function beam_second_deriv(ξ::Float64, nnodes::Int)
    if nnodes==2
        DD = Array{Float64}(4)
    else
        DD = Array{Float64}(6)
    end
    return DD
end

function elem_config_dofs(mat::AbsBeam, elem::Element)::Void
    if elem.ndim==2
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :rz, :mz)
        end
    end
    if elem.ndim==3
        for node in elem.nodes
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
            add_dof(node, :rz, :mz)
        end
    end
end

function elem_map(mat::AbsBeam, elem::Element)::Array{Int,1}
    if elem.ndim==2
        dof_keys = (:ux, :uy, :rz)
    else
        dof_keys = (:ux, :uy, :uz, :rx, :ry, :rz)
    end
    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)
end

# Return the class of element where this material can be used
client_shape_class(mat::AbsBeam) = LINE_SHAPE

function elem_dF!(mat::AbsBeam, elem::Element, dU::Array{Float64,1})
    K = elem_stiffness(mat, elem)
    return K*dU
end

function elem_and_node_vals(mat::AbsBeam, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    return node_vals, elem_vals

end
