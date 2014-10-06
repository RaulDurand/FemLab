
export PPTruss

type PPTrussIpData<:IpData
    ndim::Int
    σ::Float64
    ε::Float64
    εpa::Float64
    Δγ ::Float64

    function PPTrussIpData(ndim=3)
        this = new(ndim)
        this.σ = 0.0
        this.ε = 0.0
        this.εpa = 0.0
        this.Δγ  = 0.0
        return this
    end
end

type PPTruss<:AbsTruss
    E::Float64
    A::Float64
    σy0::Float64
    H::Float64
    new_ipdata::DataType

    function PPTruss(;E=NaN, A=NaN, sig_y=NaN, H=0.0)
        @check E>0
        @check A>0
        @check sig_y>0
        if H==0.0
            H = E*1.0e-15
            #H = 0.01
        end
        this = new(E,A,sig_y,H)
        this.new_ipdata = PPTrussIpData
        this
    end
end

function set_state(ipd::PPTrussIpData, σ=NaN, ε=NaN)
    if !isnan(σ); ipd.σ = σ end
    if !isnan(ε); ipd.ε = ε end
end

function yield_func(mat::PPTruss, ipd::PPTrussIpData, σ::Float64)
    σya = mat.σy0 + mat.H*ipd.εpa
    return abs(σ) - σya
end

function calcD(mat::PPTruss, ipd::PPTrussIpData)
    if ipd.Δγ == 0.0
        return mat.E
    else
        E, H = mat.E, mat.H
        return E*H/(E+H)
    end
end

function stress_update(mat::PPTruss, ipd::PPTrussIpData, Δε::Float64)
    E, H = mat.E, mat.H
    σini = ipd.σ
    σtr    = σini + E*Δε
    ftr    = yield_func(mat, ipd, σtr)
    ipd.Δγ = ftr>0.0? ftr/(E+H) : 0.0
    Δεp    = ipd.Δγ*sign(σtr)
    ipd.εpa += ipd.Δγ
    ipd.σ  = σtr - E*Δεp
    Δσ     = ipd.σ - σini
    ipd.ε += Δε
    return Δσ
end

function getvals(ipd::PPTrussIpData)
    return [ 
      :sa => ipd.σ,
      :ea => ipd.ε]
      #:Fa => ipd.σ*mat.A,
      #:A  => mat.A ]
end


