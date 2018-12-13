
export Truss

mutable struct TrussIpState<:IpState
    ndim::Int
    σ::Float64
    ε::Float64
    function TrussIpState(ndim::Int=3)
        this = new(ndim)
        this.σ = 0.0
        this.ε = 0.0
        return this
    end
end

mutable struct Truss<:AbsTruss
    E::Float64
    A::Float64
   xi::Float64
   ro::Float64

    function Truss(prms::Dict{Symbol,Float64})
        return  Truss(;prms...)
    end

    function Truss(;E::Number=1.0, A::Number=1.0, xi::Number=1.0, ro::Number=1.0)
        if E<=0.0; error("Invalid value for E: $E") end
        if A<=0.0; error("Invalid value for A: $A") end
        if xi<0.0; error("Invalid value for xi: $xi") end
        if ro<0.0; error("Invalid value for ro: $ro") end
        this = new(E,A,xi,ro)
        return this
    end
end

# Create a new instance of Ip data
new_ip_state(mat::Truss, ndim::Int) = TrussIpState(ndim)

function set_state(ipd::TrussIpState, σ=NaN, ε=NaN)
    if !isnan(σ); ipd.σ = σ end
    if !isnan(ε); ipd.ε = ε end
end

function stress_update(mat::Truss, ipd::TrussIpState, Δε::Float64)
    E  = mat.E
    Δσ = mat.E*Δε
    ipd.ε += Δε
    ipd.σ += Δσ
    return Δσ
end

function ip_state_vals(mat::Truss, ipd::TrussIpState)
    return Dict(
      :sa => ipd.σ,
      :ea => ipd.ε,
      :Fa => ipd.σ*mat.A,
      :A  => mat.A )
end

function calcD(mat::Truss, ips::TrussIpState)
    return mat.E
end

