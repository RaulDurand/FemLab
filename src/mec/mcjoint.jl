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
    upa::Float64  # max plastic displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size
    function MCJointIpData(ndim=3)
        this = new(ndim)
        this.σ = zeros(3)
        this.w = zeros(3)
        this.upa = 0.0
        this.Δλ  = 0.0
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
    if length(sig)==3
        ipd.σ[:] = sig
    else
        if length(sig)!=0; error("MCJointIpData: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.ε[:] = eps
    else
        if length(eps)!=0; error("MCJointIpData: Wrong size for strain array: $eps") end
    end
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
    σmax = a - upa*b
    return (σ[2]^2 + σ[3]^2)^0.5 + (σ[1]-σmax)*mat.μ
end

function yield_func2(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - upa*b
    return σ[2]^2 + σ[3]^2 - ( (σmax - σ[1])*mat.μ )^2
end

function yield_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1})
    τ = (σ[2]^2 + σ[3]^2)^0.5
    return [ mat.μ, σ[2]/τ, σ[3]/τ ]
end

function potential_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - upa*b
    if σ[1]<0.0
        return [0.0, 2*σ[2], 2*σ[3]]
    else
        τmax = σmax*mat.μ
        return [2*σ[1]/σmax^2, 2*σ[2]/τmax^2, 2*σ[3]/τmax^2]
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
        Dep  = De - De*r*v'*De/(v'*De*r - y*m*normβ(r,mat.β))
        @show Dep
        return Dep
    end
end

function stress_update(mat::MCJoint, ipd::MCJointIpData, Δw::Vect)
    #@show Δw
    σini   = ipd.σ

    # calculate De
    kn = mat.E*mat.α/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h
    De = [  kn  0.0  0.0
           0.0   ks  0.0
           0.0  0.0   ks ]

    σtr  = ipd.σ + De*Δw
    #@show σtr
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa)
    #@show Ftr

    if Ftr<0.0
        ipd.Δλ = 0.0
        ipd.σ  = σtr
    else

    maxits = 10
    upa0    = ipd.upa
    for i=1:maxits
        @show i
        # calculates Δλ
        # solving a second degree polynomial  A*Δλ^2 + B*λ + C = 0
        r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)  # r_n+1 -> r_n
        #@show r
        q = De*r
        #@show q
        a, b = calc_σmax(mat, ipd.upa)
        #@show a, b
        l  = mat.μ^2
        nr = normβ(r, mat.β)
        #@show nr
        A = q[2]^2 + q[3]^2 - b^2*nr^2*l + 2*b*nr*q[1]*l - q[1]^2*l
        #@show A
        B = -2*σtr[2]*q[2] - 2*σtr[3]*q[3] + 2*a*b*nr*l - 2*a*q[1]*l - 2*b^2*ipd.upa*nr*l + 2*b*ipd.upa*q[1]*l - 2*b*nr*σtr[1]*l + 2*σtr[1]*q[1]*l
        #@show B
        C = yield_func2(mat, ipd, σtr, ipd.upa)
        #@show C
        
        D = B^2 - 4*A*C
        #@show D
        Δλ1 = (-B + D^0.5)/(2*A)
        ipd.σ   = σtr - Δλ1*De*r
        Δupa    = Δλ1*nr
        ipd.upa = upa0 + Δupa
        #@show(Δλ1)
        #@show(ipd.σ)
        #@show(ipd.upa)
        @show yield_func(mat, ipd, ipd.σ, ipd.upa) 
        if abs( yield_func(mat, ipd, ipd.σ, ipd.upa) ) < 1e-4
            ipd.Δλ = Δλ1
            break
        else
            Δλ2 = (-B - D^0.5)/(2*A)
            ipd.σ   = σtr - Δλ2*De*r
            Δupa    = Δλ2*nr
            ipd.upa = upa0 + Δupa
            #@show(Δλ2)
            #@show(ipd.σ)
            #@show(ipd.upa)
            @show yield_func(mat, ipd, ipd.σ, ipd.upa)
            if abs( yield_func(mat, ipd, ipd.σ, ipd.upa) ) < 1e-4
                ipd.Δλ = Δλ2
                break
            #else
                #error("Δλ: something is wrong... ")
            end
        end
    end

        #ipd.σ   = σtr - ipd.Δλ*De*r
        #Δupa    = ipd.Δλ*nr
        #ipd.upa += Δupa
    end

    # update w
    ipd.w += Δw

    # calculate Δσ
    Δσ = ipd.σ - σini

    return Δσ
end

function stress_update_iter(mat::MCJoint, ipd::MCJointIpData, Δw::Vect)
    @show Δw
    σini   = ipd.σ

    # calculate De
    kn = mat.E*mat.α/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h
    De = [  kn  0.0  0.0
           0.0   ks  0.0
           0.0  0.0   ks ]

    σtr  = ipd.σ + De*Δw
    @show σtr
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa)
    @show Ftr

    if Ftr<0.0
        ipd.Δλ = 0.0
        ipd.σ  = σtr
    else
        # calculates Δλ

        maxiter = 20
        #ipd.Δλ = -1.0
        ipd.Δλ = 0.0001

        upa0 = ipd.upa

        for i=1:maxiter
            #ipd.Δλ += 0.001
            @show i
            r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)  # r_n+1 -> r_n
            #@show r
            nr = normβ(r, mat.β)
            q = De*r
            a, b = calc_σmax(mat, ipd.upa)

            num = ( (σtr[2] - ipd.Δλ*q[2])^2  + (σtr[3] - ipd.Δλ*q[3])^2 )^0.5  - (a - b*ipd.upa)*mat.μ + (σtr[1]-ipd.Δλ*q[1])*mat.μ
            den = -b*nr*mat.μ

            #num = ( (σtr[2] - ipd.Δλ*q[2])^2  + (σtr[3] - ipd.Δλ*q[3])^2 )^0.5  - a*mat.μ + b*ipd.upa*mat.μ + σtr[1]*mat.μ
            #den = mat.μ*q[1] - b*nr*mat.μ

            #@show num
            #@show den
            ipd.Δλ = num/den
            @show ipd.Δλ
            ipd.σ  = σtr - ipd.Δλ*q
            Δupa   = ipd.Δλ*nr
            ipd.upa = upa0 + Δupa
            #@show ipd.upa
            @show yield_func(mat, ipd, ipd.σ, ipd.upa)

            if abs( yield_func(mat, ipd, ipd.σ, ipd.upa) ) < 1e-4
                break
            end
        end

        #ipd.σ   = σtr - ipd.Δλ*De*r
        #Δupa    = ipd.Δλ*nr
        #ipd.upa += Δupa
    end

    # update w
    ipd.w += Δw

    # calculate Δσ
    Δσ = ipd.σ - σini

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

