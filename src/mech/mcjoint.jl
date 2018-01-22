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


mutable struct MCJointIpState<:IpState
    ndim::Int
    σ  ::Array{Float64,1}
    w  ::Array{Float64,1}
    wpa::Array{Float64,1}
    upa::Float64  # max plastic displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # element size
    function MCJointIpState(ndim=3)
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

mutable struct MCJoint<:AbsJoint
    E  ::Float64  # Young's modulus
    ν  ::Float64  # Poisson ratio
    σmax0::Float64  # tensile strength
    μ  ::Float64  # friction angle
    α  ::Float64  # elastic scale factor
    wc ::Float64  # critical openning
    ws ::Float64  # openning at inflection
    model::String # model type ("bilinear" or "hordijk")

    function MCJoint(prms::Dict{Symbol,Float64})
        return  MCJoint(;prms...)
    end

    function MCJoint(;E=NaN, nu=NaN, ft=NaN, mu=NaN, alpha=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, model="")

        if isnan(wc)
            if model == "bilinear"
                if isnan(Gf)
                    wc = round(5*GF/(1000*ft), 8)
                    ws = round(wc*0.15, 8)
                else
                    wc = round(2*(GF - Gf)/(250*ft) + 2*Gf/(1000*ft), 8)
                    ws = round(1.5*Gf/(1000*ft), 8)
                end
            elseif model == "hordijk"
                wc = round(GF/(194.7019536*ft), 8)
            end    
        end

        @assert(model=="bilinear" || model=="hordijk") 
        @assert(isnan(ws) || ws>0)
        @assert(E>0 && nu>=0 && ft>0 && mu>0 && alpha>0 && wc>0)
       
        this = new(E, nu, ft, mu, alpha, wc, ws, model)
        return this
    end
end

# Create a new instance of Ip data
new_ip_state(mat::MCJoint, ndim::Int) = MCJointIpState(ndim)

function set_state(ipd::MCJointIpState, sig=zeros(0), eps=zeros(0))
    @assert(false)
end

function calc_σmax(mat::MCJoint, upa)
    if mat.model == "bilinear"
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
        σmax = a - b*upa
    elseif mat.model == "hordijk"
        if upa < mat.wc
            z = (1 + 27*(upa/mat.wc)^3)*e^(-6.93*upa/mat.wc) - 28*(upa/mat.wc)*e^(-6.93)
        else
            z = 0
        end
        σmax = z*mat.σmax0
    end
    return σmax
end

function σmax_deriv(mat::MCJoint, upa)
    # ∂σmax/∂upa = dσmax
    if mat.model == "bilinear"
        σs = 0.25*mat.σmax0
        if upa<mat.ws
            b  = (mat.σmax0 - σs)/mat.ws
        elseif upa<mat.wc
            b  = σs/(mat.wc-mat.ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.model == "hordijk"
        if upa < mat.wc
            dz = ((81*upa^2*e^(-6.93*upa/mat.wc)/mat.wc^3) - (6.93*(1 + 27*upa^3/mat.wc^3)*e^(-6.93*upa/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0
        end
        dσmax = dz*mat.σmax0
    end
    return dσmax
end

function calc_kn_ks(ipd::MCJointIpState, mat::MCJoint)
    if mat.model == "bilinear"
        if ipd.upa<mat.ws
            α = mat.α - (mat.α - 2.0)*ipd.upa/mat.ws
        elseif ipd.upa<mat.wc
            α = 2.0 - 1.5*ipd.upa/mat.wc
        else
            α = 0.5
        end
    elseif mat.model == "hordijk"
        z = (1 + 27*(ipd.upa/mat.wc)^3)*e^(-6.93*ipd.upa/mat.wc) - 28*(ipd.upa/mat.wc)*e^(-6.93)
        if ipd.upa<mat.wc
            α = mat.α - (mat.α - 0.5)*(1-z)
        else
            α = 0.5
        end 
    end
    kn = mat.E*α/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*α/ipd.h
    return kn,ks
end


function yield_func(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1}, upa)
    σmax = calc_σmax(mat, upa)
    if ipd.ndim==3
        return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_deriv(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1})
    if ipd.ndim==3
        return [ mat.μ, σ[2]/sqrt(σ[2]^2 + σ[3]^2), σ[3]/sqrt(σ[2]^2 + σ[3]^2)]
    else
        return [ mat.μ, sign(σ[2]) ]
    end
end


function potential_derivs(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1}, upa)
    σmax = calc_σmax(mat, upa)
    if ipd.ndim==3
        if σmax>0.0
            if σ[1] >= 0.0 
                # G1:
                r = [ 2.0*σ[1]*mat.μ^2, 2.0*σ[2], 2.0*σ[3]]
            else
                # G2:
                r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
            end
        else # σmax==0.0    
                r = [ 0.0, 2*σ[2], 2*σ[3] ]
        end
    else
        if σmax>0.0
            if σ[1] >= 0.0 
                # G1:
                r = [ 2*σ[1]*mat.μ^2, 2*σ[2]]
            else
                # G2:
                r = [ 0.0, 2*σ[2] ]
            end
        else # σmax==0.0    
                r = [ 0.0, 2*σ[2] ]
        end
    end
    return r
end

function find_intersection(mat::MCJoint, ipd::MCJointIpState, F1::Float64, F2::Float64, σ0::Array{Float64,1}, Δσ::Array{Float64,1})
    @assert(F1*F2<0.0)
    α  = 0.0
    α1 = 0.0
    α2 = 1.0
    F  = F1
    maxits = 40
    for i=1:maxits
        α  = α1 + F1/(F1-F2)*(α2-α1)
        F  = yield_func(mat, ipd, σ0 + α*Δσ, ipd.upa)

        if abs(F)<1e-7 break end

        if F<0.0
            α1 = α
            F1 = F
        else
            α2 = α
            F2 = F
        end
    end

    return α
end

function calc_Δλ(mat::MCJoint, ipd::MCJointIpState, σtr::Array{Float64,1})
    μ      = mat.μ
    kn, ks = calc_kn_ks(ipd, mat)
    σmax   = calc_σmax(mat, ipd.upa)
    m      = σmax_deriv(mat, ipd.upa) 

    if σtr[1] > 0
        degree = 3
    else
        degree = 2
    end

    if degree == 2
        if ipd.ndim==3
            a = -2*ks*m*μ
            b = 2*σtr[1]*μ*ks - 2*ks*σmax*μ - m*μ
            c = (σtr[2]^2 + σtr[3]^2)^0.5 + σtr[1]*μ - σmax*μ
        else
            a = -2*m*μ*ks
            b = 2*σtr[1]*μ*ks - 2*ks*σmax*μ - m*μ
            c = abs(σtr[2]) + σtr[1]*μ - σmax*μ
        end

        Δ = b^2 - 4*a*c

        if Δ < 0
        	warn("MCJoint: Complex root x1 = $(x1) e x2 = $(x2).")
            @show(ipd.upa)
            x0 = -1
        else
        	x1 = (-b + (b^2 - 4*a*c)^0.5)/(2*a)
            x2 = (-b - (b^2 - 4*a*c)^0.5)/(2*a)
            if x1>=0 && x2>=0
            	x0 = min(x1, x2)
            else 
            	x0 = max(x1,x2)
            end
        end
    elseif degree == 3
        if ipd.ndim==3
            a = -4*kn*ks*m*(μ^3)
            b = -4*kn*ks*σmax*(μ^3) - 2*kn*m*(μ^3) - 2*ks*m*μ
            c = 2*((σtr[2]^2 + σtr[3]^2)^0.5)*kn*(μ^2) + 2*σtr[1]*μ*ks - 2*kn*σmax*(μ^3) - 2*ks*σmax*μ - m*μ
            d = (σtr[2]^2 + σtr[3]^2)^0.5 + σtr[1]*μ - σmax*μ
        else
            a = -4*kn*ks*m*(μ^3)
            b = -4*kn*ks*σmax*(μ^3) - 2*kn*m*(μ^3) - 2*ks*m*μ
            c = 2*abs(σtr[2])*kn*(μ^2) + 2*σtr[1]*μ*ks - 2*kn*σmax*(μ^3) - 2*ks*σmax*μ - m*μ
            d = abs(σtr[2]) + σtr[1]*μ - σmax*μ
        end

        f(x) = a*x^3 + b*x^2 + c*x + d
        deriv(x) = 3*a*x^2 + 2*b*x + c
        
        x0 = 1e-8
        eps = 1e-20
        i = 0
        err = eps + 1 # err > eps
        x = x0 # inicializa x

        while err > eps
            i = i + 1 # atualiza i
            x = x0 - f(x0)/deriv(x0)
            err = abs(x - x0)
            x0 = x # atualiza x0
        end

        if x0 < 0.0
            warn("MCJoint: Negative plastic multiplier Δλ = $(x0).")
            @show(ipd.upa)
        end
    end
    return x0
end

function calc_σ(mat::MCJoint, ipd::MCJointIpState, σtr::Array{Float64,1})
    μ = mat.μ
    kn, ks = calc_kn_ks(ipd, mat) 

    if ipd.ndim==3
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*ipd.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks), σtr[3]/(1 + 2*ipd.Δλ*ks)]
        end    
    else
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*ipd.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*ipd.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*ipd.Δλ*ks)]
        end    
    end
    return σ
end

function mountD(mat::MCJoint, ipd::MCJointIpState)
    kn, ks = calc_kn_ks(ipd, mat)
    σmax = calc_σmax(mat, ipd.upa)

    if ipd.ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn   0.0
                 0.0  ks  ]
    end

    if ipd.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 
        Dep  = De*0.0 
        return Dep
    else
        v    = yield_deriv(mat, ipd, ipd.σ)
        r    = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
        y    = -mat.μ # ∂F/∂σmax
        m    = σmax_deriv(mat, ipd.upa)  # ∂σmax/∂upa

        Dep  = De - De*r*v'*De/(v'*De*r - y*m*norm(r))
    	return Dep
    end
end

function stress_update(mat::MCJoint, ipd::MCJointIpState, Δw::Array{Float64,1})
    σini = copy(ipd.σ)

    μ = mat.μ
    kn, ks = calc_kn_ks(ipd, mat)
    σmax = calc_σmax(mat, ipd.upa)  

    if ipd.ndim==3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn  0.0
               0.0  ks ]
    end

    # σ trial and F trial
    σtr  = ipd.σ + De*Δw
    Ftr  = yield_func(mat, ipd, σtr, ipd.upa) 
 
    # Elastic and EP integration
    if σmax == 0.0 && ipd.w[1] >= 0.0 
        if ipd.ndim==3
            r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
            r = r1/norm(r1)
            ipd.Δλ = (σtr[1]^2- σtr[2]^2 -σtr[3]^2)/(σtr[1]*kn*r[1] - σtr[2]*ks*r[2] - σtr[3]*ks*r[3])
        else
            r1 = [ σtr[1]/kn, σtr[2]/ks ]
            r = r1/norm(r1)
            ipd.Δλ = (abs(σtr[2]) + σtr[1]*μ)/(kn*μ*r[1] + ks*r[2]*sign(σtr[2]))     
        end

        ipd.upa += ipd.Δλ
        ipd.σ = σtr - ipd.Δλ*De*r
        
        # Plastic update of w and Δσ
        ipd.w += Δw
        Δσ = ipd.σ - σini

    elseif Ftr <= 0.0
        # Pure elastic increment
        ipd.Δλ = 0.0
        ipd.σ  = copy(σtr) 
        ipd.w += Δw
        Δσ = ipd.σ - σini

    else
        # Find intersection with the yield surface
        α  = 0.0
        Fini = yield_func(mat, ipd, ipd.σ, ipd.upa)

        # Elastic increment
        if Ftr>1e-6 && Fini<0.0
            α = find_intersection(mat, ipd, Fini, Ftr, ipd.σ, De*Δw)
            ipd.w += α*Δw
            ipd.σ += α*De*Δw
        end

        σi = copy(ipd.σ)
        upai = copy(ipd.upa)

        # Elastic-plastic increment
        Δwep = (1.0-α)*Δw
        σtr = ipd.σ + De*Δwep
        ipd.Δλ = calc_Δλ(mat, ipd, σtr)
       
          
        if ipd.Δλ <0
            ipd.upa = mat.wc
            ipd.σ = zeros(ipd.ndim)
        elseif mat.model == "bilinear"
            σ = calc_σ(mat, ipd, σtr)
            r = potential_derivs(mat, ipd, σ, ipd.upa)
            upa = ipd.upa + ipd.Δλ*norm(r)

            if ipd.upa < mat.ws
                a  = mat.σmax0 
            elseif ipd.upa < mat.wc
                a  = mat.wc*0.25*mat.σmax0/(mat.wc-mat.ws)
            else
                a = 0.0
            end

            if a  == mat.σmax0 && upa < mat.ws
                ipd.σ = calc_σ(mat, ipd, σtr)
                r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
                ipd.upa += ipd.Δλ*norm(r)

            elseif a  == mat.σmax0 && upa > mat.ws   
                ipd.upa = upai
                ipd.σ = σi
                n = 30
                Δwinc = Δwep/n

                for i= 1:n
                    σtr = ipd.σ + De*Δwinc
                    ipd.Δλ = calc_Δλ(mat, ipd, σtr)   
                    ipd.σ = calc_σ(mat, ipd, σtr)
                    r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
                    ipd.upa += ipd.Δλ*norm(r)
                end
            else
                ipd.σ = calc_σ(mat, ipd, σtr)
                r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
                ipd.upa += ipd.Δλ*norm(r)
            end    
        elseif mat.model == "hordijk"
            ipd.σ = calc_σ(mat, ipd, σtr)
            r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
            ipd.upa += ipd.Δλ*norm(r)
        end

        # Return to surface:
         F  = yield_func(mat, ipd, ipd.σ, ipd.upa)

        #if false 
        if F>1e-2
            ipd.upa = upai
            ipd.σ = σi
            n = 100
            Δwinc = Δwep/n

            for i= 1:n
            	σtr = ipd.σ + De*Δwinc
                ipd.Δλ = calc_Δλ(mat, ipd, σtr)
                ipd.σ = calc_σ(mat, ipd, σtr)
                r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)
                ipd.upa += ipd.Δλ*norm(r)
            end
        end
        
        #  Plastic update of w and Δσ
        ipd.w += Δwep
        Δσ = ipd.σ - σini
    end

    return Δσ
end

function ip_state_vals(mat::MCJoint, ipd::MCJointIpState)
    if ipd.ndim == 3
       return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] ,
          :upa   => ipd.upa
          )
    else
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :upa   => ipd.upa
          )
    end
end
