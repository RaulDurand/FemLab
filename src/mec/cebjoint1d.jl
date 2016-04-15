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


export CEBJoint1D

type CEBJoint1DIpData<:IpData
    ndim::Int
    σ  ::Array{Float64,1}
    ε  ::Array{Float64,1}
    τy ::Float64
    unload::Bool
    σc::Float64
    function CEBJoint1DIpData(ndim=3)
        this = new(ndim)
        this.σ = zeros(3)
        this.ε = zeros(3)
        this.τy     = 0.0
        this.unload = true
        return this
    end
end

type CEBJoint1D<:AbsJoint1D
    h ::Float64
    ks::Float64
    kn::Float64
    s1::Float64
    s2::Float64
    s3::Float64
    τres::Float64

    new_ipdata::DataType

    function CEBJoint1D(prms::Dict{Symbol,Float64})
        return  CEBJoint1D(;prms...)
    end

    function CEBJoint1D(;ks=NaN, kn=NaN, TauR=NaN, s1=NaN, s2=NaN, s3=NaN, h=NaN, A=NaN, dm=NaN)
        @assert ks>=0
        @assert ks*s1>TauR
        @assert s2>s1
        @assert s3>s2

        if isnan(kn) 
            kn = ks
        end
        @assert kn>=0

        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^0.5
            else
                h = pi*dm
            end
        end
        @assert h>0

        this = new(h, ks, kn, s1, s2, s3, TauR)
        this.new_ipdata = CEBJoint1DIpData
        return this
    end
end

function Tau(mat::CEBJoint1D, s::Float64)
    ss = abs(s)
    #@assert s>=0
    if ss<mat.s1
        return mat.ks*s
    elseif ss<mat.s2
        return mat.ks*mat.s1*sign(s)
    elseif ss<mat.s3
        τmax = mat.ks*mat.s1
        #return (mat.τres-τmax)/(mat.s3-mat.s2)*(mat.s3-ss)*sign(s)
        return (mat.τres + (mat.τres-τmax)/(mat.s3-mat.s2)*(ss-mat.s3))*sign(s)
        
    else
        return mat.τres*sign(s)
    end
end

function deriv(mat::CEBJoint1D, ipd::CEBJoint1DIpData, s::Float64)
    ss = abs(s)
    
    #check for unloading... !
    #τ = ipd.σ[1]
    #f = yield_func(mat, ipd, τ, ss)


    if ss<mat.s1
        return mat.ks
    elseif ss<mat.s2
        return 0.0
    elseif ss<mat.s3
        τmax = mat.ks*mat.s1
        return (mat.τres-τmax)/(mat.s3-mat.s2)
    else
        return 0.0
    end
end

function set_state(ipd::CEBJoint1DIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.σ[:] = sig
    else
        if length(sig)!=0; error("CEBJoint1DIpData: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.ε[:] = eps
    else
        if length(eps)!=0; error("CEBJoint1DIpData: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::CEBJoint1D, ipd::CEBJoint1DIpData)

    # inital value for ipd.τy
    if ipd.τy==0.0
        ipd.τy = mat.ks*mat.s1
    end

    s  = abs(ipd.ε[1])
    kh = deriv(mat, ipd, s)

    ks = ipd.unload? mat.ks : mat.ks*kh/(mat.ks + kh)

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

function yield_func(mat::CEBJoint1D, ipd::CEBJoint1DIpData, τ::Float64, s::Float64)
    return abs(τ) - ipd.τy
end

function stress_update(mat::CEBJoint1D, ipd::CEBJoint1DIpData, Δε::Vect)
    ks = mat.ks
    kn = mat.kn
    Δs = Δε[1]      # relative displacement
    τini = ipd.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial
    s    = ipd.ε[1]

    ftr  = yield_func(mat, ipd, τtr, abs(s+Δs) )

    if ftr<0.0
        τ = τtr 
        ipd.unload = true
        ftr  = yield_func(mat, ipd, τtr, abs(s+Δs) )
    else
        τ   = Tau(mat, s+Δs)
        Δτ  = τ - τini
        Δse  = Δτ/ks
        Δsp  = abs(Δs - Δse)
        ipd.τy = min(ipd.τy, abs(Tau(mat, s+Δs)))
        ipd.unload = abs(ipd.ε[1]+Δs) - abs(s) >= 0.0
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

function getvals(mat::CEBJoint1D, ipd::CEBJoint1DIpData)
    return Dict(
      :ur   => ipd.ε[1] ,
      :tau  => ipd.σ[1] ,
      )
end

