export DruckerPrager
export set_state

type DruckerPragerIpData<:IpData
    ndim::Int
    σ::Array{Float64,1}
    ε::Array{Float64,1}
    Δγ::Float64
    pl::Bool
    function DruckerPragerIpData(ndim=3) 
        this = new(ndim)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this
    end
end

type DruckerPrager<:Mechanical
    E::Float64
    ν::Float64
    α::Float64
    κ::Float64
    new_ipdata::DataType

    function DruckerPrager(;E=NaN, nu=0.0, alpha=0.0, kappa=0.0)
        @check E>=0.0
        @check 0<=nu<0.5
        @check alpha>=0.0
        @check kappa>=0.0

        this = new(E, nu, alpha, kappa)
        this.new_ipdata = DruckerPragerIpData

        return this
    end
end

function set_state(ipd::DruckerPragerIpData; sig=zeros(0), eps=zeros(0))
    sq2 = √2.
    mdl = [1., 1., 1., sq2, sq2, sq2]
    if length(sig)==6
        ipd.σ[:] = sig.*mdl
    else
        if length(sig)!=0; error("DruckerPrager: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*mdl
    else
        if length(eps)!=0; error("DruckerPrager: Wrong size for strain array: $eps") end
    end
end

function yield_func(mat::DruckerPrager, σ::Tensor2)
    j1  = J1(σ)
    j2d = J2D(σ)
    p   = j1/3
    α,κ = mat.α, mat.κ
    return α*j1 + √j2d - κ
end

function yield_deriv(mat::DruckerPrager, σ::Tensor2) # df/dσ
    j1  = J1(σ)
    if j1>mat.T
        return tI
    end

    j2d = J2D(σ)
    if j2d==0.0
        return tI
    end

    α,κ = mat.α, mat.κ
    s   = dev(mat.σ)
    return α*tI + 0.5*s/√j2d  # TODO check
end


function mountD(mat::DruckerPrager, ipd::DruckerPragerIpData)
    E, ν = mat.E, mat.ν
    c = E/((1+ν)*(1-2*ν))
    [ c*(1.-ν)      c*ν        c*ν            0.0            0.0           0.0
          c*ν   c*(1.-ν)       c*ν            0.0            0.0           0.0
          c*ν       c*ν    c*(1.-ν)           0.0            0.0           0.0
          0.0       0.0        0.0    c*(1.-2.*ν)            0.0           0.0
          0.0       0.0        0.0            0.0    c*(1.-2.*ν)           0.0
          0.0       0.0        0.0            0.0            0.0    c*(1.-2.*ν) ]
end

function mountDep(mat::DruckerPrager, ipd::DruckerPragerIpData)
    De = mountD(mat, ipd)
    if !ipd.plas
        return De
    end

    v = yield_deriv(ipd.σ)

    return De - (De∷v) ⊗ (v∷De) / (v∷De∷v)

end

function stress_update(mat::DruckerPrager, ipd::DruckerPragerIpData, deps::Array{Float64,1})
    D    = mount_D(mat, ipd)
    dsig = D*deps
    ipd.ε += deps
    ipd.σ += dsig
    dsig
end

function getvals(mat::DruckerPrager, ipd::DruckerPragerIpData)
    σ  = ipd.σ
    ε  = ipd.ε
    ndim = ipd.ndim
    sr2  = 2.0^0.5

    if ndim==2;
        return [
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :sxy => σ[4]/sr2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :exy => ε[4]/sr2,
          :s_m => sum(σ[1:3])/3.0 ]
      else
        return [
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :sxy => σ[4]/sr2,
          :syz => σ[5]/sr2,
          :sxz => σ[6]/sr2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :exy => ε[4]/sr2,
          :eyz => ε[5]/sr2,
          :exz => ε[6]/sr2,
          :s_m => sum(σ[1:3])/3.0 ]
      end
end

