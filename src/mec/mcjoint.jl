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
        #b = 1.0
        b  = -σs/(mat.wc-mat.ws)/1e6
        #b = -1e-1
    end
    return a, b
end

#function yield_func(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    #a, b = calc_σmax(mat, upa)
    #σmax = a - b*upa
    #return abs(σ[2]) + abs(σ[3]) + (σ[1]-σmax)*mat.μ
#end

#function yield_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1})
    #return [ mat.μ, sign(σ[2]), sign(σ[3]) ]
#end

function yield_func(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - b*upa
    return (σ[2]^2 + σ[3]^2)^0.5 + (σ[1]-σmax)*mat.μ
end

function yield_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1})
    τ = (σ[2]^2 + σ[3]^2)^0.5
    return [ mat.μ, σ[2]/τ, σ[3]/τ ]
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
        Dep  = De - De*r*v'*De/(v'*De*r - y*m*dot(l,r))

        return Dep
    end
end

function find_intersection(mat::MCJoint, ipd::MCJointIpData, F1::Float64, F2::Float64, σ0::Array{Float64,1}, Δσ::Array{Float64,1})
    @assert(F1*F2<0.0)

    α  = 0.0
    α1 = 0.0
    α2 = 1.0
    F  = F1
    while abs(F)>1e-5
        α  = α1 + F1/(F1-F2)*(α2-α1)
        F  = yield_func(mat, ipd, σ0 + α*Δσ, ipd.upa)
        if F<0.0
            α1 = α
            F1 = F
        else
            α2 = α
            F2 = F
        end
        #@show F
    end

    return α
end

    using PyPlot
function stress_update(mat::MCJoint, ipd::MCJointIpData, Δw::Vect)
    σini   = copy(ipd.σ)

    # calculate De
    kn = mat.E*mat.α/ipd.h

    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h
    De = [  kn  0.0  0.0
           0.0   ks  0.0
           0.0  0.0   ks ]

    σtr  = ipd.σ + De*Δw
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa)


    if Ftr<=0.0
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr)

        # update w
        ipd.w += Δw

        # calculate Δσ
        Δσ = ipd.σ - σini
    else
        # Pure elastic increment
        α  = 0.0
        # Find intersection with the yield surface
        Fini = yield_func(mat, ipd, ipd.σ, ipd.upa)
        if Ftr>1e-6 && Fini<0.0
            α = find_intersection(mat, ipd, Fini, Ftr, ipd.σ, De*Δw)
            # Update w and Δσ up to the intersection
            ipd.w += α*Δw
            ipd.σ += α*De*Δw
        end

        # Elastic-plastic increment
        Δwep = (1.0-α)*Δw

        y   = -mat.μ # dF/dσmax
        upa = ipd.upa
        r   = zeros(3)

        F  = yield_func(mat, ipd, ipd.σ, ipd.upa)

        nincs = 10
        Δwi = Δwep/nincs
        for i=1:nincs
            σtr = ipd.σ + De*Δwi
            v   = yield_derivs(mat, ipd, ipd.σ)
            r   = potential_derivs(mat, ipd, ipd.σ, ipd.upa)

            if ipd.Δλ==0.0
                l = zeros(3)
            else
                l = ipd.wpa .* [1.0, mat.β, mat.β ] / ipd.upa # dupa/dwpa
            end

            a, b = calc_σmax(mat, upa)
            m    = -b   # dσmax/dupa

            ipd.Δλ = dot(v, De*Δwi)/( dot(v,De*r) - y*m*dot(l,r) )

            #q = De*r
            #ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/ ( q[2]*sign(σtr[2]) + q[3]*sign(σtr[3]) +q[1]*mat.μ - b*dot(l,r)*mat.μ  )

            if ipd.Δλ<0.0
                warn("MCJoint: Negative plastic multiplier Δλ = $(ipd.Δλ)")
            end

            Δwpa = ipd.Δλ*r
            ipd.wpa += Δwpa
            ipd.upa  = normβ(ipd.wpa, mat.β)

            ipd.σ = σtr - ipd.Δλ*De*r

            # Return to surface:
            F  = yield_func(mat, ipd, ipd.σ, ipd.upa)

            #if false 
            if F>1e-6
                Δσ = ipd.σ
                # The return algorithm needs improving
                σ0 = [-0.000001, 0.0, 0.0]
                F0 = yield_func(mat, ipd, σ0, ipd.upa)
                @assert(F0*F<0.0)
                α  = find_intersection(mat, ipd, F0, F, σ0, Δσ)

                ipd.σ = σ0 + α*Δσ
            end

        end

        # Elastic plastic update of w and Δσ
        ipd.w += Δwep
        Δσ = ipd.σ - σini
        F  = yield_func(mat, ipd, ipd.σ, ipd.upa)
    end

    #@show Δσ

    return Δσ
end


function getvals(mat::MCJoint, ipd::MCJointIpData)
    return Dict(
      :w1  => ipd.w[1] ,
      :w2  => ipd.w[2] ,
      :w3  => ipd.w[3] ,
      :s1  => ipd.σ[1] ,
      :s2  => ipd.σ[2] ,
      :s3  => ipd.σ[3] ,
      #:tau  => ipd.σ[1] ,
      #:sigc => ipd.σc   ,
      #:taumax => τmax   ,
      :upa   => ipd.upa
      )
end

