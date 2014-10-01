export set_mat

abstract Mechanical<:Material

function config_dofs(::Mechanical, elem::Element)
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        if elem.ndim==3; add_dof(node, :uz, :fz) end
    end
end

function get_map(elem::Element)
    get_map(elem.mat, elem)
end

function get_map(::Mechanical, elem::Element)
    dof_keys = (:ux, :uy, :uz)[1:elem.ndim]
    vcat([ [node.dofdict[key].eq_id for key in dof_keys] for node in elem.nodes]...)
end

function set_facet_bc(mat::Mechanical, oelem::Element, face::Face, key::Symbol, val::Float64)
    fnodes = face.nodes
    fshape = face.shape
    ndim   = oelem.ndim

    if key in (:fx, :fy, :fz)
        error("Boundary condition $key is not allowed in a face; consider tx, ty, tz or tn")
    end

    if !(key in (:ux, :uy, :uz, :tx, :ty, :tz, :tn))
        error("Boundary condition $key is not applicable in a face of an element with material $(typeof(mat))")
    end

    if key in (:tz, :uz) && ndim==2
        error("Boundary condition $key is not applicable in a 2D analysis")
    end

    # Apply the boundary conditions
    if key in (:ux, :uy, :uz)
        for node in fnodes
            dof = node.dofdict[key]
            dof.prescU = true
            dof.bryU   = val
        end
        return
    end

    # Force boundary condition
    nfnodes = length(fnodes)

    # Calculate the face coordinates matrix
    C = getcoords(fnodes, ndim)

    # Calculate the vector with values to apply
    V = zeros(ndim)
    if key == :tx ; V[1] = val end
    if key == :ty ; V[2] = val end
    if key == :tz ; V[3] = val end

    # Calculate the nodal values
    F = zeros(nfnodes, ndim)

    f_ips = get_ip_coords(fshape)

    for i=1:size(f_ips,1)
        R = vec(f_ips[i,:])
        w = R[end]
        N = shape_func(fshape, R)
        D = deriv_func(fshape, R)
        J = D*C
        nJ = norm2(J)

        if key == :tn && ndim==2
            n = [J[1,2], -J[1,1]]
            V = val*n/norm(n)
        end
        if key == :tn && ndim==3
            n = cross(J[1,:], J[2,:])
            V = val*n/norm(n)
        end

        F += N*V'*(nJ*w)
    end
    #@show F

    # Setting bc into nodes
    for (i,node) in enumerate(fnodes)
        node.dofdict[:fx].bryF += F[i,1]
        node.dofdict[:fy].bryF += F[i,2]
        if ndim==3; node.dofdict[:fz].bryF += F[i,3] end
    end
end

function mount_B(::Mechanical, elem::Element, R::Vect, C::Matx, B::Matx)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    D = deriv_func(elem.shape, R)
    J = D*C
    dNdX = inv(J)*D
    detJ = det(J)
    sqr2 = 2.0^0.5
    B[:] = 0.0
    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[1,i]
            B[2,2+j*ndim] = dNdX[2,i]
            B[4,1+j*ndim] = dNdX[2,i]/sqr2; B[4,2+j*ndim] = dNdX[1,i]/sqr2
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[1,i]
            dNdy = dNdX[2,i]
            dNdz = dNdX[3,i]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,1+j*ndim] = dNdy/sqr2;   B[4,2+j*ndim] = dNdx/sqr2
            B[5,2+j*ndim] = dNdz/sqr2;   B[5,3+j*ndim] = dNdy/sqr2
            B[6,1+j*ndim] = dNdz/sqr2;   B[6,3+j*ndim] = dNdx/sqr2
        end
    end

    return detJ
end


function stiff(elem::Element)
    stiff(elem.mat, elem)
end

function stiff(::Mechanical, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    mat    = elem.mat
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)

    for ip in elem.ips
        detJ = mount_B(mat, elem, ip.R, C, B)
        D    = mount_D(mat, ip.data)
        coef = detJ*ip.w
        K   += B'*D*B*coef
    end
    return K
end

function update!(elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    update!(elem.mat, elem, DU, DF)
end

function update!(::Mechanical, elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    map    = get_map(elem)
    mat    = elem.mat

    dF = zeros(nnodes*ndim)
    dU = DU[map]
    B  = zeros(6, nnodes*ndim)

    C = getcoords(elem)
    for ip in elem.ips
        detJ = mount_B(mat, elem, ip.R, C, B)
        deps = B*dU
        dsig = stress_update(mat, ip.data, deps)
        coef = detJ*ip.w
        dF += B'*dsig*coef
    end

    # Update global vector
    DF[map] += dF

end


function getvals(::Mechanical, ips::Array{Ip,1})
    all_ip_vals = [ get_vals(ip.data) for ip in ips ]
    labels      = keys(all_ip_vals[1])
    nips        = length(ips) 

    # matrix with all ip values (nip x nvals)
    IP = vcat([ [values(all_ip_vals[i])...] for i=1:nips]...)
    avg_vals = mean(IP,1)

    # filling elem values dict
    Dict(labels, avg_vals)

end

function node_and_elem_vals(mat::Mechanical, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    # Elem vals
    all_ip_vals = [ getvals(mat, ip.data) for ip in elem.ips ]
    labels      = keys(all_ip_vals[1])
    nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    E = extrapolator(elem.shape, nips)
    N = E*IP # (nnodes x nvals)


    # Filling nodal and elem vals
    for (i,key) in enumerate(labels)
        node_vals[key] = N[:,i]
        #elem_vals[key] = mean(IP[:,i])
    end

    return node_vals, elem_vals

end


include("solid.jl")
include("truss.jl")
include("pptruss.jl")
include("joint1d.jl")
include("mcjoint1d.jl")

