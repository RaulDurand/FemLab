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


export MCJoint1D

type MCJoint1DIpData<:IpData
    ndim::Int
    σ  ::Array{Float64,1}
    ε  ::Array{Float64,1}
    ωpa::Float64
    Δγ :: Float64
    σc :: Float64
    function MCJoint1DIpData(ndim=3)
        this = new(ndim)
        this.σ = zeros(3)
        this.ε = zeros(3)
        this.ωpa = 0.0
        this.Δγ  = 0.0
        this.σc  = 0.0
        return this
    end
end

type MCJoint1D<:AbsJoint1D
    ks::Float64
    kn::Float64
    h ::Float64
    c ::Float64
    μ ::Float64
    kh::Float64
    new_ipdata::DataType

    function MCJoint1D(;ks=NaN, kn=NaN, h=NaN, A=NaN, dm=NaN, c=NaN, C=NaN, mu=NaN, phi=NaN, kh=0.0)
        @check ks>=0
        @check (h>=0 || A>0 || dm>0)
        @check (c>=0 || C>=0)
        @check (phi>=0 || mu>=0)

        if isnan(kn)
            kn = ks
        end
        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^2.0
            else
                h = pi*dm
            end
        end
        @check h>=0
        if isnan(c)
            c = C
        end
        if isnan(mu)
            mu = tan(phi) # μ = tan(φ)
        end

        this = new(ks, kn, h, c, mu, kh)
        this.new_ipdata = MCJoint1DIpData
        return this
    end
end

function set_state(ipd::MCJoint1DIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.σ[:] = sig
    else
        if length(sig)!=0; error("MCJoint1DIpData: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.ε[:] = eps
    else
        if length(eps)!=0; error("MCJoint1DIpData: Wrong size for strain array: $eps") end
    end
end

function mountD(mat::MCJoint1D, ipd::MCJoint1DIpData)
    ks = ipd.Δγ==0.0? mat.ks : mat.ks*mat.kh/(mat.ks + mat.kh)
    kn = mat.kn
    if ipd.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end

function yield_func(mat::MCJoint1D, ipd::MCJoint1DIpData, τ::Float64)
    σc = ipd.σc>0? 0.0 : abs(ipd.σc)
    c  = mat.c
    kh = mat.kh
    μ  = mat.μ
    return f = abs(τ) - (c + kh*ipd.ωpa + μ*σc)
end

function stress_update(mat::MCJoint1D, ipd::MCJoint1DIpData, Δε::Vect)
    ks = mat.ks
    kn = mat.kn
    kh = mat.kh
    Δω = Δε[1] # relative displacement
    τini = ipd.σ[1]
    τtr  = τini + ks*Δω
    ftr  = yield_func(mat, ipd, τtr)

    if ftr<0.0
        ipd.Δγ = 0.0
        τ = τtr 
    else
        ipd.Δγ = ftr/(ks+kh)
        Δωp    = ipd.Δγ*sign(τtr)
        τ      = τtr - ks*Δωp # correcting first term
        ipd.ωpa += ipd.Δγ
    end

    # update ε
    ipd.ε += Δε

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δε
    Δσ[1] = Δτ

    # update Δσ
    ipd.σ += Δσ

    return Δσ
end


function getvals(mat::MCJoint1D, ipd::MCJoint1DIpData)
    τmax = mat.c + abs(ipd.σc)*mat.μ
    return [ 
      :ur   => ipd.ε[1] ,
      :tau  => ipd.σ[1] ,
      :sigc => ipd.σc   ,
      :taumax => τmax   ,
      :w_pa   => ipd.ωpa
      ]
end

