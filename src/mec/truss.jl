
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
    return Dict(
      :sa => ipd.σ,
      :ea => ipd.ε,
      :Fa => ipd.σ*mat.A,
      :A  => mat.A )
end

function calcD(mat::Truss, ips::TrussIpData)
    return mat.E
end

function elem_jacobian(mat::AbsTruss, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    A = mat.A
    E = mat.E
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

function update!(mat::AbsTruss, elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    map    = get_map(elem)
    A      = mat.A

    dF = zeros(nnodes*ndim)
    dU = DU[map]
    C = getcoords(elem)
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
