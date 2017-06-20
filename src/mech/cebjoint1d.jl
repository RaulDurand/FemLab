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


mutable struct CEBJoint1DIpState<:IpState
    ndim::Int
    σ  ::Array{Float64,1}
    ε  ::Array{Float64,1}
    τy ::Float64      # max stress
    sy ::Float64      # accumulated relative displacement
    elastic::Bool
    function CEBJoint1DIpState(ndim=3)
        this = new(ndim)
        this.σ = zeros(3)
        this.ε = zeros(3)
        this.τy = 0.0
        this.sy = 0.0
        this.elastic = false
        return this
    end
end

mutable struct CEBJoint1D<:AbsJoint1D
    τmax:: Float64
    τres:: Float64
    s1  :: Float64
    s2  :: Float64
    s3  :: Float64
    α   :: Float64
    β   :: Float64
    ks  :: Float64
    kn  :: Float64
    h   :: Float64

    function CEBJoint1D(prms::Dict{Symbol,Float64})
        return  CEBJoint1D(;prms...)
    end

    function CEBJoint1D(;TauM=NaN, TauR=NaN, s1=NaN, s2=NaN, s3=NaN, alpha=NaN, beta=NaN, kn=NaN, ks=NaN, h=NaN, A=NaN, dm=NaN)
        @assert s1>0
        @assert s2>s1
        @assert s3>s2
        @assert ks>0

        # Estimate the perimeter h
        if isnan(h)
            if A>0
                h = 2.0*√(A*pi)
            else
                @assert dm>0
                h = pi*dm
            end
        end
        @assert h>0

        # Estimate TauM if not provided
        if isnan(TauM)
            TauM = ks*s1
        end
        @assert TauM>TauR

        # Define alpha if not provided
        if isnan(alpha); alpha = 1.0 end
        @assert 0.0<=alpha<=1.0

        # Define beta if not provided
        if isnan(beta); beta= 1.0 end
        @assert beta>=0.0

        # Estimate kn if not provided
        if isnan(kn)
            kn = ks
        end
        @assert kn>0

        this = new(TauM, TauR, s1, s2, s3, alpha, beta, ks, kn, h)
        return this
    end
end

# Create a new instance of Ip data
new_ip_state(mat::CEBJoint1D, ndim::Int) = CEBJoint1DIpState(ndim)

function Tau(mat::CEBJoint1D, sy::Float64)
    if sy<mat.s1
        return mat.τmax*(sy/mat.s1)^mat.α
    elseif sy<mat.s2
        return mat.τmax
    elseif sy<mat.s3
        return mat.τmax - (mat.τmax-mat.τres)*((sy-mat.s2)/(mat.s3-mat.s2))^mat.β
    else
        return mat.τres
    end
end

function deriv(mat::CEBJoint1D, ipd::CEBJoint1DIpState, sy::Float64)
    if sy==0.0
        const S1_FACTOR = 0.01
        sy = S1_FACTOR*mat.s1   # to avoid undefined derivative
    end

    if sy<=mat.s1
        return mat.τmax/mat.s1*(sy/mat.s1)^(mat.α-1)
    elseif sy<mat.s2
        return mat.ks/1000
        #return 0.0
    elseif sy<mat.s3
        return -(mat.τmax-mat.τres)/(mat.s3-mat.s2)*((sy-mat.s2)/(mat.s3-mat.s2))^(mat.β-1)
    else
        #return 0.0
        return mat.ks/1000
    end
end

function set_state(ipd::CEBJoint1DIpState, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.σ[:] = sig
    else
        if length(sig)!=0; error("CEBJoint1DIpState: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.ε[:] = eps
    else
        if length(eps)!=0; error("CEBJoint1DIpState: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::CEBJoint1D, ipd::CEBJoint1DIpState)
    if ipd.elastic
        ks = mat.ks
    else
        ks = deriv(mat, ipd, ipd.sy)
    end

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

function yield_func(mat::CEBJoint1D, ipd::CEBJoint1DIpState, τ::Float64)
    return abs(τ) - ipd.τy
end

function stress_update(mat::CEBJoint1D, ipd::CEBJoint1DIpState, Δε::Vect)
    ks = mat.ks
    kn = mat.kn
    Δs = Δε[1]      # relative displacement
    τini = ipd.σ[1] # initial shear stress
    τtr  = τini + ks*Δs # elastic trial

    ftr  = yield_func(mat, ipd, τtr)

    if ftr<0.0
        τ = τtr 
        ipd.elastic = true
    else
        Δsy     = (abs(τtr)-ipd.τy)/mat.ks
        ipd.sy += Δsy
        ipd.τy  = Tau(mat, ipd.sy)
        τ  = ipd.τy*sign(τtr)
        Δτ = τ - τini
        ipd.elastic = false
    end

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δε
    Δσ[1] = Δτ

    # update ε and σ
    ipd.ε += Δε
    ipd.σ += Δσ

    return Δσ
end

function ip_state_vals(mat::CEBJoint1D, ipd::CEBJoint1DIpState)
    return Dict(
      :ur   => ipd.ε[1] ,
      :tau  => ipd.σ[1] ,
      )
end

