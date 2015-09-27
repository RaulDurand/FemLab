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

export DPCap
export set_state

type DPCapIpData<:IpData
    ndim::Int
    σ::Tensor2
    ε::Tensor2
    εpa::Float64
    Δγ::Float64
    function DPCapIpData(ndim=3) 
        this = new(ndim)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.εpa = 0.0
        this.Δγ  = 0.0
        this
    end
end

type DPCap<:Mechanical
    E::Float64
    ν::Float64
    α::Float64
    κ::Float64
    ρ::Float64
    H::Float64
    De::Tensor4
    new_ipdata::DataType

    function DPCap(;E=NaN, nu=0.0, alpha=0.0, kappa=0.0, rho=0.0, H=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert alpha>=0.0
        @assert kappa>0.0
        @assert rho>0.0
        @assert H>=0.0

        this     = new(E, nu, alpha, kappa, rho, H)
        this.De  = mount_De(E,nu) # elastic tensor
        this.new_ipdata = DPCapIpData
        return this 
    end
end

function set_state(ipd::DPCapIpData; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig.*V2M
    else
        if length(sig)!=0; error("DPCap: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*V2M
    else
        if length(eps)!=0; error("DPCap: Wrong size for strain array: $eps") end
    end
end

function yield_func(mat::DPCap, ipd::DPCapIpData, σ::Tensor2)
    j1  = J1(σ)
    j2d = J2D(σ)
    α,κ = mat.α, mat.κ
    ρ   = mat.ρ
    H   = mat.H
    εpa = ipd.εpa
    return √j2d - (κ-α*j1)*(1. - j1^2./(ρ+H*εpa)^2.)
end

function mount_D(mat::DPCap, ipd::DPCapIpData)
    α   = mat.α
    H   = mat.H

    De = mat.De
    if ipd.Δγ==0.0
        return De
    end

    j2d = J2D(ipd.σ)
    if j2d != 0.0
        s  = dev(ipd.σ) 
        su = s/norm(s)
        α, H  = mat.α, mat.H
        κ     = mat.κ
        ρ  = mat.ρ
        j1 = J1(ipd.σ)
        r  = α*(1. - j1^2./(ρ+H*ipd.εpa)^2.) - 2.*(κ-α*j1)*j1/(ρ+H*ipd.εpa)^2.  #df/dj1
        V  = r*tI + su/√2 # df/dσ
        N  = V
        Nu = N/norm(N)
    else # apex
        j1 = J1(ipd.σ)
        if j1>0.0
            Nu = 1./√3.*tI
        else
            Nu = -1./√3.*tI
        end
        V  = Nu
    end

    return De - inner(De,Nu) ⊗ inner(V,De) / (inner(V,De,Nu) + H)

end

function secant_root(f::Function, x0::Float64)
    eps = 1.0e-8
    err = eps + 1.0
    x   = x0
    h   = 1.0
    x0 = 0.001
    #@show f(x)
    while abs(f(x))>eps
        dfdx = (f(x+h)-f(x-h))/(2.0*h)
        x    = x0 - f(x0)/dfdx
        h    = x - x0
        x0   = x
        #@show x
        #@show f(x)
    end 
        #@show x
        #@show f(x)
    #x = 0.0
    #for i=1:100
        #@show x
        #@show(f(x))
        #x += 0.002
    #end
    return x
end


function stress_update(mat::DPCap, ipd::DPCapIpData, Δε::Array{Float64,1})
    σini   = ipd.σ
    σtr    = ipd.σ + inner(mat.De, Δε)
    ftr    = yield_func(mat, ipd, σtr)


    if ftr <= 0.0
        # elastic
        ipd.Δγ = 0.0
        ipd.σ  = σtr
    else
        # plastic 
        K, G  = mat.E/(3.*(1.-2.*mat.ν)), mat.E/(2.*(1.+mat.ν))
        α, H  = mat.α, mat.H
        κ     = mat.κ
        ρ     = mat.ρ
        j1    = J1(ipd.σ) # warning: evaluated at step n and not at n+1
        #@show  j1
        r     = α*(1. - j1^2./(ρ+H*ipd.εpa)^2.) - 2.*(κ-α*j1)*j1/(ρ+H*ipd.εpa)^2.  #df/dj1
        n     = 1./√(3.*r*r+0.5)
        j1tr  = J1(σtr)
        #@show  j1tr
        j2dtr = J2D(σtr)

        if j1 < -mat.ρ# return to max cap
            #j1 = -mat.rho
            fn1 = (Δγ) -> begin 
                j1n1 = j1tr - 3*√3*Δγ*K;
                (κ-α*j1n1) * (1. - j1n1^2./(ρ+H*ipd.εpa+Δγ*H)^2.)
            end
            ipd.Δγ    = secant_root(fn1, 0.0)
            j1    = j1tr -3*√3*ipd.Δγ*K;
            ipd.σ = j1/3.*tI
        else
            if √j2dtr - ipd.Δγ*n*G > 0.0  # conventional return
                # find Δγ
                fn1 = (Δγ) -> begin 
                    j1n1 = j1tr-9.*Δγ*r*n*K;
                    √j2dtr - Δγ*n*G - (κ-α*j1n1) * (1. - j1n1^2./(ρ+H*ipd.εpa+Δγ*H)^2.)
                end
                #@show j1
                ipd.Δγ    = secant_root(fn1, 0.0)
                jj1 = j1
                j1    = j1tr - 9.*ipd.Δγ*r*n*K
                #if ipd.Δγ<0.0
                if false
                    @show jj1
                    @show ipd.Δγ
                    @show √j2dtr
                    @show n
                    @show G
                    @show √j2dtr - ipd.Δγ*n*G # conventional return
                    #g=-0.01
                    #for i=1:100
                        #g += 0.001 
                        #println(g, "  ", fn1(g))
                    #end
                    exit()
                end
                m     = 1. - ipd.Δγ*n*G/√j2dtr
                ipd.σ = m*dev(σtr) + j1/3.*tI
            else # return to apex
                #exit()
                #kkk=1
                sgn = j1tr>0.0? -1.0: 1.0
                # find Δγ
                fn1 = (Δγ) -> begin 
                    j1n1 = j1tr + sgn*3*√3*Δγ*K;
                    (κ-α*j1n1) * (1. - j1n1^2./(ρ+H*ipd.εpa+Δγ*H)^2.)
                end
                ipd.Δγ    = secant_root(fn1, 0.0)
                #@show ipd.Δγ
                j1    = j1tr + sgn*3*√3*ipd.Δγ*K;
                #@show j1
                ipd.σ = j1/3.*tI
            end
        end

        ipd.εpa += ipd.Δγ
    end

    ipd.ε += Δε
    Δσ     = ipd.σ - σini
    return Δσ
end

function getvals(mat::DPCap, ipd::DPCapIpData)
    σ  = ipd.σ
    ε  = ipd.ε
    ndim = ipd.ndim
    sr2  = √2.

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
          :dg  => ipd.Δγ,
          :j1  => trace(σ),
          :srj2d => √J2D(σ),
          :p   => trace(σ)/3. ]
      end
end

