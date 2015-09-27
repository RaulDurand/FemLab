
export Truss

abstract AbsTruss<:Mechanical

type TrussIpData<:IpData
    ndim::Int
    σ::Float64
    ε::Float64
    function TrussIpData(ndim::Int=3)
        this = new(ndim)
        this.σ = 0.0
        this.ε = 0.0
        return this
    end
end

type Truss<:AbsTruss
    E::Float64
    A::Float64
    new_ipdata::DataType

    function Truss(prms::Dict{Symbol,Float64})
        return  Truss(;prms...)
    end

    function Truss(;E::Float64=1.0, A::Float64=1.0)
        if E<=0.0; error("Invalid value for E: $E") end
        if A<=0.0; error("Invalid value for A: $A") end
        this = new(E,A)
        this.new_ipdata = TrussIpData
        this
    end
end

function set_state(ipd::TrussIpData, σ=NaN, ε=NaN)
    if !isnan(σ); ipd.σ = σ end
    if !isnan(ε); ipd.ε = ε end
end

function stress_update(mat::Truss, ipd::TrussIpData, Δε::Float64)
    E  = mat.E
    Δσ = mat.E*Δε
    ipd.ε += Δε
    ipd.σ += Δσ
    return Δσ
end

function getvals(mat::Truss, ipd::TrussIpData)
    return [ 
      :sa => ipd.σ,
      :ea => ipd.ε,
      :Fa => ipd.σ*mat.A,
      :A  => mat.A ]
end


function mount_B(::AbsTruss, elem::Element, R::Vect, C::Matx, B::Matx)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    dNdR = deriv_func(elem.shape, R)
    J = dNdR*C
    detJ = norm(J)
    B[:] = 0.0

    for i in 1:nnodes
        for j=1:ndim
            B[1,j+(i-1)*ndim] = dNdR[1,i]*J[j]/detJ^2.0
        end
    end

    return detJ
end

function calcD(mat::Truss, ips::TrussIpData)
    return mat.E
end

function elem_jacobian(mat::AbsTruss, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    A = elem.mat.A
    C = getcoords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)

    for ip in elem.ips
        detJ = mount_B(elem.mat, elem, ip.R, C, B)
        D    = calcD(mat,ip.data)

        coef = D*A*detJ*ip.w
        K   += B'*B*coef
    end
    return K
end

function update!(mat::AbsTruss, elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    map    = get_map(elem)
    A      = mat.A

    dF = zeros(nnodes*ndim)
    dU = DU[map]
    B  = zeros(1, nnodes*ndim)

    C = getcoords(elem)
    for ip in elem.ips
        detJ = mount_B(elem.mat, elem, ip.R, C, B)
        deps = (B*dU)[1]
        dsig = stress_update(mat, ip.data, deps)
        coef = A*detJ*ip.w
        dF  += B'*dsig*coef
    end

    # Update global vector
    DF[map] += dF
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
