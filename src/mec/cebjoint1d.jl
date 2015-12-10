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
    ωpa::Float64
    Δγ ::Float64
    σc ::Float64 #not used
    function CEBJoint1DIpData(ndim=3)
        this = new(ndim)
        this.σ = zeros(3)
        this.ε = zeros(3)
        this.ωpa = 0.0
        this.Δγ  = 0.0
        return this
    end
end

type CEBJoint1D<:AbsJoint1D
    Tau::Array{Float64,2}
    h ::Float64
    ks::Float64
    kn::Float64

    new_ipdata::DataType

    function CEBJoint1D(prms::Dict{Symbol,Float64})
        return  CEBJoint1D(;prms...)
    end

    #function CEBJoint1D(;τmax=NaN, τres=NaN, s1=NaN, s2=NaN, s3=NaN, kn=NaN, h=NaN, kn=Nan, A=NaN, dm=NaN)
    function CEBJoint1D(;tau=[], h=NaN, kn=NaN, A=NaN, dm=NaN)
        @assert length(tau)>0
        @assert kn>=0

        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^0.5
            else
                h = pi*dm
            end
        end
        @assert h>0

        ks = (tau[2,2] - tau[1,2])/(tau[2,1] - tau[1,1])

        this = new(tau, h, ks, kn)
        this.new_ipdata = CEBJoint1DIpData
        return this
    end
end

function interpolate(M::Array{Float64,2}, val::Float64)
    #@show val
    #@show M
    @assert val>=M[1,1] && val<=M[end,1]
    for i=2:size(M,1)
        if val<=M[i,1]
            return M[i-1,2] + (M[i,2]-M[i-1,2])/(M[i,1]-M[i-1,1])*(val - M[i-1,1])
        end
    end
    #@show "Whatt"
    return NaN
end

function deriv(M::Array{Float64,2}, val::Float64)
    @assert val>=M[1,1] && val<=M[end,1]
    for i=2:size(M,1)
        if val<=M[i,1]
            return (M[i,2]-M[i-1,2])/(M[i,1]-M[i-1,1])
        end
    end
    #@show "Whatt2"
    return NaN
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

function mountD(mat::CEBJoint1D, ipd::CEBJoint1DIpData)
    ω  = abs(ipd.ε[1])
    kh = deriv(mat.Tau, ω)
    #@show ω
    #@show kh
    ks = ipd.Δγ==0.0? mat.ks : mat.ks*kh/(mat.ks + kh)
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

function yield_func(mat::CEBJoint1D, ipd::CEBJoint1DIpData, τ::Float64, ω::Float64)
    #ω = abs(ipd.ε[1])
    #@show ω
    #@show interpolate(mat.Tau, ω) 
    return f = abs(τ) - interpolate(mat.Tau, ω) 
end

function bisection(f::Function, a, b)
    if f(a) == 0.0
        return a
    end

    while f(a)*f(b) >= 0
        b *= 2.0
    end

    eps = 1e-8
    n = floor(Int, log(2, (b-a)/eps) + 1)
    x = 0.0

    for i=1:n
        x = (a+b)/2
        if f(a)*f(x) < 0
            b = x
        else
            a = x
        end
    end

    return x
end

function stress_update(mat::CEBJoint1D, ipd::CEBJoint1DIpData, Δε::Vect)
    ks = mat.ks
    kn = mat.kn
    Δω = Δε[1] # relative displacement
    τini = ipd.σ[1]
    τtr  = τini + ks*Δω
    ftr  = yield_func(mat, ipd, τtr, ipd.ωpa+abs(Δω) )
    #@show τtr
    #@show ftr

    if ftr<0.0
        ipd.Δγ = 0.0
        τ = τtr 
    else
        f(Δγ) = abs(τtr) - ks*Δγ - interpolate(mat.Tau, ipd.ωpa +Δγ)
        ipd.Δγ  = bisection(f, 0.0, 1.0)
        #@show Δγ
        #@show f(Δγ)
        Δωp = ipd.Δγ*sign(τtr)
        τ   = τtr - ks*Δωp # correcting first term
        ipd.ωpa += ipd.Δγ
        #@show ipd.ωpa
        #@show yield_func(mat, ipd, τ, ipd.ωpa )
    end
    #@show τ

    # update ε
    #@show Δε
    ipd.ε += Δε

    # calculate Δσ
    Δτ = τ - τini
    Δσ = kn*Δε
    Δσ[1] = Δτ

    # update Δσ
    ipd.σ += Δσ

    #@show ipd.ε[1]
    #@show ipd.σ[1]

    return Δσ
end

function getvals(mat::CEBJoint1D, ipd::CEBJoint1DIpData)
    τmax = mat.c + abs(ipd.σc)*mat.μ
    return [ 
      :ur   => ipd.ε[1] ,
      :tau  => ipd.σ[1] ,
      :sigc => ipd.σc   ,
      :taumax => τmax   ,
      :w_pa   => ipd.ωpa
      ]
end

