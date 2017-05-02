
export Truss

abstract AbsTruss<:Mechanical

# Return the class of element where this material can be used
client_elem_class(mat::AbsTruss) = :LINE

function elem_jacobian(mat::AbsTruss, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    A = mat.A
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)
    J  = Array(Float64, 1, ndim)

    for ip in elem.ips
        dNdR = deriv_func(elem.shape, ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)

        # mount B
        B[:] = 0.0
        for i in 1:nnodes
            for j=1:ndim
                B[1,j+(i-1)*ndim] = dNdR[1,i]*J[j]/detJ^2.0
            end
        end

        E    = calcD(mat,ip.data)
        coef = E*A*detJ*ip.w
        @gemm K += coef*B'*B
    end
    return K
end

function update!(mat::AbsTruss, elem::Element, dU::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    A      = mat.A

    dF = zeros(nnodes*ndim)
    C  = getcoords(elem)
    B  = zeros(1, nnodes*ndim)
    J  = Array(Float64, 1, ndim)
    for ip in elem.ips
        dNdR = deriv_func(elem.shape, ip.R)
        @gemm J = dNdR*C
        detJ = norm(J)

        # mount B
        B[:] = 0.0
        for i in 1:nnodes
            for j=1:ndim
                B[1,j+(i-1)*ndim] = dNdR[1,i]*J[j]/detJ^2.0
            end
        end

        deps = (B*dU)[1]
        dsig = stress_update(mat, ip.data, deps)
        coef = A*detJ*ip.w
        dF  += coef*B'*dsig
    end

    return dF

end

function node_and_elem_vals(mat::AbsTruss, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    # Elem vals
    all_ip_vals = [ getvals(mat, ip.data) for ip in elem.ips ]
    # completing with axial forces
    for ip_val in all_ip_vals
        ip_val[:A ] = mat.A
        ip_val[:Fa] = mat.A*ip_val[:sa]
    end
    labels      = keys(all_ip_vals[1])
    nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    E = extrapolator(elem.shape, nips)
    N = E*IP # (nnodes x nvals)


    # Filling nodal and elem vals
    for (i,key) in enumerate(labels)
        #node_vals[key] = N[:,i] # No nodal values for truss
        elem_vals[key] = mean(IP[:,i])
    end

    return node_vals, elem_vals

end
