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
        this.σ = zeros(ndim)
        this.w = zeros(ndim)
        this.wpa = zeros(ndim)
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
        b  = -σs/(mat.wc-mat.ws)/1.0e8 # important in pull tests
    end
    return a, b
end


function yield_func(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - b*upa
    if ipd.ndim==3
        return √(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1})
    if ipd.ndim==3
        τ = √(σ[2]^2 + σ[3]^2)
        if τ==0.0
            return [ 1., 0., 0.]
        else
            return [ mat.μ, σ[2]/τ, σ[3]/τ ]
        end
    else
        τ = σ[2]
        if τ==0.0
            return [ 1., 0]
        else
            return [ mat.μ, sign(τ) ]
        end
    end
end


function potential_derivs(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - upa*b
    if ipd.ndim==3
        if σ[1]>0.0
            # G1:
            σm = √(σ[2]^2 + σ[3]^2 + σmax^2*mat.μ^2)
            if σm == 0.0 return [ 1., 0., 0. ] end
            r  = [ σ[1]*mat.μ^2/σm, σ[2]/σm, σ[3]/σm ]
            return r/norm(r)
        else
            # G2:
            τ = √(σ[2]^2 + σ[3]^2)
            if τ != 0.0
                r = [ 0.0, σ[2]/τ, σ[3]/τ ]
                r/norm(r)
            else
                r = [ 1., 0., 0.]
            end
            return r
        end
    else
        if σ[1]>0.0
            # G1:
            σm = √(σ[2]^2 + σmax^2*mat.μ^2)
            if σm == 0.0 return [ 1., 0. ] end
            r = [ σ[1]*mat.μ^2/σm, σ[2]/σm ]
            return r/norm(r)
        else
            # G2:
            τ = σ[2]
            if τ != 0.0
                r = [ 0., sign(τ) ]
            else
                r = [ 1., 0. ]
            end
            return r
        end

    end
end


function mountD(mat::MCJoint, ipd::MCJointIpData)
    kn = mat.E*mat.α/ipd.h

    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h
    if ipd.ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn   0.0
                 0.0  ks  ]
    end

    if ipd.Δλ==0.0 # Elastic 
        return De
    else
        v    = yield_derivs(mat, ipd, ipd.σ)
        r    = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
        y    = -mat.μ # dF/dσmax
        a, b = calc_σmax(mat, ipd.upa)
        m    = -b   # dσmax/dupa

        Dep  = De - De*r*v'*De/(v'*De*r - y*m*norm(r))

        return Dep
    end
end


function find_intersection(mat::MCJoint, ipd::MCJointIpData, F1::Float64, F2::Float64, σ0::Array{Float64,1}, Δσ::Array{Float64,1})
    @assert(F1*F2<0.0)

    α  = 0.0
    α1 = 0.0
    α2 = 1.0
    F  = F1
    while abs(F)>1e-7
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


function SecantRoot(f::Function, x0::Float64, x1::Float64, eps::Float64=1e-6)
    # initial values
    x  = x0
    f0 = f(x0)
    maxits = 20
    i=0
    for i=1:maxits
        f1 = f(x1)
        d  = (f1-f0)/(x1-x0) 
        x  = x1 - f1/d
        err = abs(x-x1)

        if err < eps
            return x, true
        end

        # updating
        x0 = x1
        x1 = x
        f0 = f1
    end

    return x, false
end


function stress_update(mat::MCJoint, ipd::MCJointIpData, Δw::Array{Float64,1})
    σini = copy(ipd.σ)
    μ    = mat.μ

    # calculate De
    kn = mat.E*mat.α/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h

    if ipd.ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn  0.0
               0.0  ks ]
    end

    # Trial
    σtr  = ipd.σ + De*Δw
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa)


    # Elastic and EP integration
    if Ftr<=0.0
        # Pure elastic increment
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr)

        # update w
        ipd.w += Δw

        # calculate Δσ
        Δσ = ipd.σ - σini
    else
        # Elastic-plastic increment
        upa = ipd.upa

        nincs = 1
        Δwi = Δw/nincs
        for i=1:nincs
            σtr = ipd.σ + De*Δwi
            v   = yield_derivs(mat, ipd, ipd.σ)
            r   = potential_derivs(mat, ipd, ipd.σ, ipd.upa)

            a, b = calc_σmax(mat, upa)
            m    = -b   # dσmax/dupa

            if ipd.ndim==3
                ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ^2-b*μ)
            else
                ipd.Δλ = yield_func(mat, ipd, σtr, ipd.upa)/(kn*μ^2 - b*μ + ks*r[2]*sign(σtr[2]))
            end


            if ipd.Δλ<0.0
                warn("MCJoint: Negative plastic multiplier Δλ = $(ipd.Δλ)")
                @show i
                @show yield_func(mat, ipd, σtr, ipd.upa) # if <0 then elastic
                @show (kn*mat.μ^2-b*mat.μ)
            end

            Δwpa = ipd.Δλ*r
            ipd.upa += norm(Δwpa)

            ipd.σ = σtr - ipd.Δλ*De*r

            # Return to surface:
            F  = yield_func(mat, ipd, ipd.σ, ipd.upa)

            #if false 
            if F>1e-6
                # The return algorithm needs improving
                Δσ = ipd.σ
                #σ0 = [-0.000001, 0.0, 0.0]
                σ0 = [-0.000001, 0.0, 0.0][1:ipd.ndim]
                F0 = yield_func(mat, ipd, σ0, ipd.upa)
                @assert(F0*F<0.0)
                α  = find_intersection(mat, ipd, F0, F, σ0, Δσ)

                ipd.σ = σ0 + α*Δσ
            end

        end

        # Elastic plastic update of w and Δσ
        ipd.w += Δw
        Δσ = ipd.σ - σini
        F  = yield_func(mat, ipd, ipd.σ, ipd.upa)
    end

    return Δσ
end



function getvals(mat::MCJoint, ipd::MCJointIpData)
    if ipd.ndim == 2
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :upa   => ipd.upa
          )
    else
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] ,
          :upa   => ipd.upa
          )
    end
end

