##############################################################################
#    FemLab - Finite Element Library                                         #
#    Copyright (C) 2014 Raul Durand <raul.durand at gmail.com>               #
#                                                                            #
#    This file is part of FemLab.                                            #
#                                                                            #
#    FemLab is free software: you can redistribute it and/or modify          #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    any later version.                                                      #
#                                                                            #
#    FemLab is distributed in the hope that it will be useful,               #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with FemLab.  If not, see <http://www.gnu.org/licenses/>.         #
##############################################################################

export DPConc
export set_state

type DPConcIpData<:IpData
    ndim::Int64
    σ::Tensor2
    ε::Tensor2
    εpa::Float64
    Δγ::Float64
    function DPConcIpData(ndim=3) 
        this = new(ndim)
        this.σ   = zeros(6)
        this.ε   = zeros(6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end

type DPConc<:Mechanical
    ν::Float64
    α::Float64
    κ::Float64
    H::Float64
    fc::Float64
    εc::Float64
    #De::Tensor4
    new_ipdata::DataType

    function DPConc(prms::Dict{Symbol,Float64})
        return DPConc(;prms...)
    end

    function DPConc(;nu=0.0, alpha=0.0, kappa=0.0, H=0.0, fc=0.0, eps_c=0.0)
        @check 0.0<=nu<0.5
        @check alpha>=0.0
        @check kappa>0.0
        @check H>=0.0

        this     = new(nu, alpha, kappa, H, fc, eps_c)
        this.new_ipdata = DPConcIpData
        return this 
    end
end

function nlE(fc::Float64, εc::Float64, ε::Array{Float64,1})
    #ε = sum(ε[1:3])
    ε = maxabs(ε[1:3])
    #εt = 0.0005

    # traction
    #if ε>εt
        #return 1000.0
    #end

    # compression
    if abs(ε)>=εc
        return 1.0
    end
    return 2*fc*(εc-ε)/εc^2
end

function set_state(ipd::DPConcIpData; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig.*V2M
    else
        if length(sig)!=0; error("DPConc: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*V2M
    else
        if length(eps)!=0; error("DPConc: Wrong size for strain array: $eps") end
    end
end

function yield_func(mat::DPConc, ipd::DPConcIpData, σ::Tensor2)
    j1  = J1(σ)
    j2d = J2D(σ)
    α,κ = mat.α, mat.κ
    H   = mat.H
    εpa = ipd.εpa
    return α*j1 + √j2d - κ - H*εpa
end

function mount_D(mat::DPConc, ipd::DPConcIpData)
    α   = mat.α
    H   = mat.H

    E  = nlE(mat.fc, mat.εc, ipd.ε)
    De = mount_De(E, mat.ν) # elastic tensor

    if ipd.Δγ==0.0
        return De
    end

    j2d = J2D(ipd.σ)
    if j2d != 0.0
        s  = dev(ipd.σ) 
        su = s/norm(s)
        V  = α*tI + su/√2 # df/dσ
        N  = V
        Nu = N/norm(N)
    else # apex
        Nu = 1./√3.*tI
        V  = Nu
    end

    return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)

end

function stress_update(mat::DPConc, ipd::DPConcIpData, Δε::Array{Float64,1})
    E  = nlE(mat.fc, mat.εc, ipd.ε)
    De = mount_De(E, mat.ν) # elastic tensor

    σini   = ipd.σ
    σtr    = ipd.σ + inner(De, Δε)
    ftr    = yield_func(mat, ipd, σtr)

    if ftr < 1.e-8
        # elastic
        ipd.Δγ = 0.0
        ipd.σ  = σtr
    else
        # plastic 
        K, G  = E/(3.*(1.-2.*mat.ν)), E/(2.*(1.+mat.ν))
        #K, G  = mat.E/(3.*(1.-2.*mat.ν)), mat.E/(2.*(1.+mat.ν))
        α, H  = mat.α, mat.H
        n     = 1./√(3.*α*α+0.5)
        j1tr  = J1(σtr)
        j2dtr = J2D(σtr)

        if √j2dtr - ipd.Δγ*n*G > 0.0 # conventional return
            ipd.Δγ = ftr/(9*α*α*n*K + n*G + H)
            j1     = j1tr - 9*ipd.Δγ*α*n*K
            m      = 1. - ipd.Δγ*n*G/√j2dtr
            ipd.σ  = m*dev(σtr) + j1/3.*tI
        else # return to apex
            κ      = mat.κ
            ipd.Δγ = (α*j1tr-κ-H*ipd.εpa)/(3*√3*α*K + H)
            j1     = j1tr - 3*√3*ipd.Δγ*K
            ipd.σ  = j1/3.*tI
        end

        ipd.εpa += ipd.Δγ

    end

    ipd.ε += Δε
    Δσ     = ipd.σ - σini
    return Δσ
end

function getvals(mat::DPConc, ipd::DPConcIpData)
    σ  = ipd.σ
    ε  = ipd.ε
    ndim = ipd.ndim
    j1   = trace(σ)
    sr2  = √2.
    srj2d = √J2D(σ)
    #pl_r  = srj2d/(mat.κ- mat.α*j1)
    #pl_r  = srj2d/j1

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
          :p   => sum(σ[1:3])/3.0 ]
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
          :ev  => trace(ε),
          :epa => trace(ipd.εpa),
          :dg  => ipd.Δγ,
          :j1  => j1,
          :srj2d => srj2d,
          :p   => trace(σ)/3.0
          #:pl_r=> pl_r
          ]
      end
end