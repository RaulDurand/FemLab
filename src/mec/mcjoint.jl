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


export MCJoint

type MCJointIpData<:IpData
    ndim::Int
    σ  ::Array{Float64,1}
    w  ::Array{Float64,1}
    wpa::Array{Float64,1}
    upa::Float64  # max plastic displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size
    function MCJointIpData(ndim=3)
        this = new(ndim)
        this.σ = zeros(3)
        this.w = zeros(3)
        this.wpa = zeros(3)
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h  = 0.0
        return this
    end
end

type MCJoint<:AbsJoint
    E  ::Float64  # Young's modulus
    ν  ::Float64  # Poisson ratio
    σmax0::Float64  # tensile strength
    μ  ::Float64  # friction angle
    α  ::Float64  # elastic scale factor
    β  ::Float64  # coupling factor for openning and sliding
    wc ::Float64  # critical openning
    ws ::Float64  # openning at inflection
    new_ipdata::DataType

    function MCJoint(prms::Dict{Symbol,Float64})
        return  MCJoint(;prms...)
    end

    function MCJoint(;E=NaN, nu=NaN, ft=NaN, mu=NaN, alfa=NaN, beta=NaN, wc=NaN, ws=NaN)
        this = new(E, nu, ft, mu, alfa, beta, wc, ws)
        this.new_ipdata = MCJointIpData
        return this
    end
end

function set_state(ipd::MCJointIpData, sig=zeros(0), eps=zeros(0))
    @assert(false)
end

function calc_σmax(mat::MCJoint, upa)
    σs = 0.25*mat.σmax0
    if upa<mat.ws
        a  = mat.σmax0 
        b  = (mat.σmax0 - σs)/mat.ws
    elseif upa<mat.wc
        a  = mat.wc*σs/(mat.wc-mat.ws)
        b  = σs/(mat.wc-mat.ws)
    else
        a = 0.0
        b = 0.0
    end
    return a, b
end

function yield_func(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - b*upa
    return abs(σ[2]) + abs(σ[3]) + (σ[1]-σmax)*mat.μ
end

function yield_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1})
    return [ mat.μ, sign(σ[2]), sign(σ[3]) ]
end

function potential_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - upa*b
    if σ[1]<=0.0 || σmax == 0.0
        r = [0.0, 2*σ[2], 2*σ[3]]
        return r/norm(r)
    else
        τmax = σmax*mat.μ
        r = [2*σ[1]/σmax^2, 2*σ[2]/τmax^2, 2*σ[3]/τmax^2]
        return r/norm(r)
    end
end

function normβ(V::Array{Float64,1}, β::Float64)
    return (V[1]^2 + β*V[2]^2 + β*V[3]^2)^0.5
end

function mountD(mat::MCJoint, ipd::MCJointIpData)
    kn = mat.E*mat.α/ipd.h

    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h
    De = [  kn  0.0  0.0
           0.0   ks  0.0
           0.0  0.0   ks ]

    #@show De

    if ipd.Δλ==0.0 # Elastic 
        return De
    else
        v    = yield_derivs(mat, ipd, ipd.σ)
        r    = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
        y    = -mat.μ # dF/dσmax
        a, b = calc_σmax(mat, ipd.upa)
        m    = -b   # dσmax/dupa

        l    = ipd.wpa .* [1.0, mat.β, mat.β ] / ipd.upa # dupa/dwpa
        #@show l
        Dep  = De - De*r*v'*De/(v'*De*r - y*m*dot(l,r))

        return Dep
    end
end

    using PyPlot
function stress_update(mat::MCJoint, ipd::MCJointIpData, Δw::Vect)
    #@show Δw
    #@show ipd.w[1]
    σini   = ipd.σ

    # calculate De
    kn = mat.E*mat.α/ipd.h

    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h
    De = [  kn  0.0  0.0
           0.0   ks  0.0
           0.0  0.0   ks ]

    σtr  = ipd.σ + De*Δw
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa)


    if Ftr<0.0
        ipd.Δλ = 0.0
        ipd.σ  = σtr

        # update w
        ipd.w += Δw

        # calculate Δσ
        Δσ = ipd.σ - σini
    else
        # Pure elastic increment

        F1 = yield_func(mat, ipd, ipd.σ, ipd.upa)
        α  = 0.0
        # Find intersection with the yield surface
        if abs(F1)>1e-6
            #F1 = -abs(F1)
            F0 = F1
            F2 = Ftr
            F = 1.0
            i = 0

            @assert(F1*F2<0.0)
            while F>1e-5
                i = i+1
                α  = F1/(F1 - F2)
                F  = yield_func(mat, ipd, ipd.σ + α*De*Δw, ipd.upa)
                if F<0.0
                    F1 = F
                else
                    F2 = F
                end
            end

            # Update w and Δσ up to the intersection
            ipd.w += α*Δw
            ipd.σ += α*De*Δw
        end

        # Elastic-plastic increment
        Δwep = (1.0-α)*Δw

        #ipd.Δλ = 0.0
        σtr = ipd.σ + De*Δwep
        y   = -mat.μ # dF/dσmax
        σ   = ipd.σ
        upa = ipd.upa
        wpa = ipd.wpa
        r   = zeros(3)
        ipd.Δλ = 0.0001

        maxit = 40
        for i=1:maxit
            v    = yield_derivs(mat, ipd, σ)
            r    = potential_derivs(mat, ipd, σ, upa)
            Δwpa = ipd.Δλ*r
            wpa  = ipd.wpa + Δwpa
            upa  = normβ(wpa, mat.β)


            a, b = calc_σmax(mat, upa)
            m    = -b   # dσmax/dupa
            l    = wpa .* [1.0, mat.β, mat.β ] / upa # dupa/dwpa


            #@show 
            @show v
            @show l
            #@show upa
            #@show ipd.upa + ipd.Δλ*dot(l,r)
            #@show v
            @show a,b
            @show r

            ipd.Δλ = dot(v, De*Δwep)/( dot(v,De*r) - y*m*dot(l,r) )
            σ      = σtr - ipd.Δλ*De*r

            @show ipd.Δλ
            F  = yield_func(mat, ipd, σ, upa)
            #@show F
            if abs(F)<1e-10
                break
            end
        end
        @assert(i<maxit)

        Δwpa = ipd.Δλ*r
        ipd.wpa += Δwpa
        ipd.upa  = normβ(ipd.wpa, mat.β)

        ipd.σ = σtr - ipd.Δλ*De*r

        # Elastic plastic update of w and Δσ
        ipd.w += Δwep
        Δσ = ipd.σ - σini
        #@show Δwpa
        #@show ipd.upa
        F  = yield_func(mat, ipd, ipd.σ, ipd.upa)
        #@show F

        #exit()
    end


    return Δσ
end


function getvals(mat::MCJoint, ipd::MCJointIpData)
    return Dict(
      :w1  => ipd.w[1] ,
      :w2  => ipd.w[2] ,
      :w3  => ipd.w[3] ,
      :σ1  => ipd.σ[1] ,
      :σ2  => ipd.σ[2] ,
      :σ3  => ipd.σ[3] ,
      #:tau  => ipd.σ[1] ,
      #:sigc => ipd.σc   ,
      #:taumax => τmax   ,
      :upa   => ipd.upa
      )
end

