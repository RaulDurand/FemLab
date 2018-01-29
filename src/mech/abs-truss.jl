
abstract type AbsTruss<:Mechanical end

# Return the class of element where this material can be used
client_shape_class(mat::AbsTruss) = LINE_SHAPE

function elem_stiffness(mat::AbsTruss, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    A = mat.A
    C = elem_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(1, ndim)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
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

function elem_mass(mat::AbsTruss, elem::Element) 
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    mat    = elem.mat
    ro = elem.mat.ro
    A = mat.A
    C = elem_coords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    J  = Array{Float64}(1, ndim)
    N = Array{Float64}(ndim, ndim*nnodes)
    

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        Ni = elem.shape.func(ip.R) #encontrei em shape.jl FemMesh
        setNt(ndim,Ni,N)    

        @gemm J = dNdR*C
        detJ = norm(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ro*A*detJ*ip.w
        @gemm M += coef*N'*N #falta multiplicar pela espessura.. perguntar onde esta.. acho que ja esta inlcuida pois precisa K

    end
    return M
end


function setNt(ndim::Int,Ni::Vect, N::Matx)
    nnodes = length(Ni) 
    N[:] = 0.0

    if ndim==2
        for i in 1:nnodes
            j = i-1
            N[1,1+j*ndim] = Ni[i]
            N[2,2+j*ndim] = Ni[i]
        end
    elseif ndim==3
        for i in 1:nnodes
            j    = i-1
            N[1,1+j*ndim] = Ni[i]
            N[2,2+j*ndim] = Ni[i]
            N[3,3+j*ndim] = Ni[i]
       end
    else
        for i in 1:nodes
            j = i-1
            N[1,1+j*ndim] = Ni[i]
        end    
    end
    

end


function elem_dF!(mat::AbsTruss, elem::Element, dU::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    A      = mat.A

    dF = zeros(nnodes*ndim)
    C  = elem_coords(elem)
    B  = zeros(1, nnodes*ndim)
    J  = Array{Float64}(1, ndim)
    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
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

function elem_and_node_vals(mat::AbsTruss, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    for key in (:vx, :vy, :vz)[1:ndim]
        node_vals[key] = [node.dofdict[Symbol("u"string(key)[2:2])].V for node in elem.nodes]
    end

    for key in (:ax, :ay, :az)[1:ndim]
        node_vals[key] = [node.dofdict[Symbol("u"string(key)[2:2])].A for node in elem.nodes]
    end

    # Elem vals
    all_ip_vals = [ ip_state_vals(mat, ip.data) for ip in elem.ips ]
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

function apply_elem_bc(mat::AbsTruss, elem::Element, key::Symbol, val::Float64)
    
    ndim   = elem.ndim
    shape  = elem.shape
    nodes = elem.nodes
    A = mat.A

    if !(key in (:gx, :gy, :gz))
        error("Boundary condition $key is not allowed in an elem; consider gx, gy or gz")
    end
   
    nnodes = length(nodes)

    # Calculate the face coordinates matrix
    C = nodes_coords(nodes, ndim)

    # Calculate the vector with values to apply
    V = zeros(ndim)
    if key == :gx ; V[1] = val end
    if key == :gy ; V[2] = val end
    if key == :gz ; V[3] = val end

    # Calculate the nodal values
    F = zeros(nnodes, ndim)

    ips = get_ip_coords(shape)

    for i=1:size(ips,1)
        R = vec(ips[i,:])
        w = R[end]
        N = shape.func(R)
        D = shape.deriv(R)
        J = D*C
        nJ = norm2(J)

        F += A*N*V'*(nJ*w)
    end

    # Setting bc into nodes
    for (i,node) in enumerate(elem.nodes)
        node.dofdict[:fx].bryF =  Expr(:call, :+, node.dofdict[:fx].bryF, F[i,1])
        node.dofdict[:fy].bryF =  Expr(:call, :+, node.dofdict[:fy].bryF, F[i,2])
        if ndim==3; node.dofdict[:fz].bryF = Expr(:call, :+, node.dofdict[:fz].bryF, F[i,3]) end
    end

end

