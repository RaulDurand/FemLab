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
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    upa::Float64  # effective plastic relative displacement
    Δλ ::Float64  # plastic multiplier
    h  ::Float64  # characteristic length from bulk elements
    function MCJointIpState(ndim=3)
        this = new(ndim)
        this.σ = zeros(ndim)
        this.w = zeros(ndim)
        this.upa = 0.0
        this.Δλ  = 0.0
        this.h  = 0.0
        return this
    end
end

mutable struct MCJoint<:AbsJoint
    E  ::Float64  # Young's modulus
    ν  ::Float64  # Poisson ratio
    σmax0::Float64  # tensile strength (internal variable)
    μ  ::Float64  # tangent of friction angle
    α  ::Float64  # maximum elastic relative displacements
    wc ::Float64  # critical crack opening
    ws ::Float64  # openning at inflection (where the curve slope changes)
    softcurve::String # softening curve model ("linear" or bilinear" or "hordijk")

    function MCJoint(prms::Dict{Symbol,Float64})
        return  MCJoint(;prms...)
    end

    function MCJoint(;E=NaN, nu=NaN, ft=NaN, mu=NaN, alpha=NaN, wc=NaN, ws=NaN, GF=NaN, Gf=NaN, softcurve="bilinear")

        if isnan(wc)
        	if softcurve == "linear"
        		 wc = round(2*GF/(1000*ft), 10)
            elseif softcurve == "bilinear"
                if isnan(Gf)
                    wc = round(5*GF/(1000*ft), 10)
                    ws = round(wc*0.15, 10)
                else
                    wc = round((8*GF- 6*Gf)/(1000*ft), 10)
                    ws = round(1.5*Gf/(1000*ft), 10)
                end
            elseif softcurve == "hordijk"
                wc = round(GF/(194.7019536*ft), 10)
            end    
        end


        @assert(softcurve=="linear" || softcurve=="bilinear" || softcurve=="hordijk") 
        @assert(isnan(ws) || ws>0)
        @assert(E > 0 && nu >= 0 && ft > 0 && mu > 0 && alpha > 0 && wc > 0)

        @show(ws)
        @show(wc)
        this = new(E, nu, ft, mu, alpha, wc, ws, softcurve)
        return this
    end
end

# Create a new instance of Ip data
new_ip_state(mat::MCJoint, ndim::Int) = MCJointIpState(ndim)

function set_state(ipd::MCJointIpState, sig=zeros(0), eps=zeros(0))
    @assert(false)
end

function yield_func(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1})
    σmax = calc_σmax(mat, ipd, ipd.upa)
    if ipd.ndim == 3
        return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_deriv(mat::MCJoint, ipd::MCJointIpState)
    if ipd.ndim == 3
        return [ mat.μ, ipd.σ[2]/sqrt(ipd.σ[2]^2 + ipd.σ[3]^2), ipd.σ[3]/sqrt(ipd.σ[2]^2 + ipd.σ[3]^2)]
    else
        return [ mat.μ, sign(ipd.σ[2]) ]
    end
end

function potential_derivs(mat::MCJoint, ipd::MCJointIpState, σ::Array{Float64,1})
    if ipd.ndim == 3
            if σ[1] >= 0.0 
                # G1:
                r = [ 2.0*σ[1]*mat.μ^2, 2.0*σ[2], 2.0*σ[3]]
            else
                # G2:
                r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
            end
    else
            if σ[1] >= 0.0 
                # G1:
                r = [ 2*σ[1]*mat.μ^2, 2*σ[2]]
            else
                # G2:
                r = [ 0.0, 2*σ[2] ]
            end
    end
    return r
end

function calc_σmax(mat::MCJoint, ipd::MCJointIpState, upa::Float64)
	if mat.softcurve == "linear"
		if upa < mat.wc
            a = mat.σmax0
            b = mat.σmax0/mat.wc
		else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*upa
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.σmax0
        if upa < mat.ws
            a  = mat.σmax0 
            b  = (mat.σmax0 - σs)/mat.ws
        elseif upa < mat.wc
            a  = mat.wc*σs/(mat.wc-mat.ws)
            b  = σs/(mat.wc-mat.ws)
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*upa
    elseif mat.softcurve == "hordijk"
        if upa < mat.wc
            z = (1 + 27*(upa/mat.wc)^3)*e^(-6.93*upa/mat.wc) - 28*(upa/mat.wc)*e^(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.σmax0
    end
    return σmax
end

function σmax_deriv(mat::MCJoint, ipd::MCJointIpState, upa::Float64)
    # ∂σmax/∂upa = dσmax
    if mat.softcurve == "linear"
		if upa < mat.wc
            b = mat.σmax0/mat.wc
		else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "bilinear"
        σs = 0.25*mat.σmax0
        if upa < mat.ws
            b  = (mat.σmax0 - σs)/mat.ws
        elseif upa < mat.wc
            b  = σs/(mat.wc-mat.ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softcurve == "hordijk"
        if upa < mat.wc
            dz = ((81*upa^2*e^(-6.93*upa/mat.wc)/mat.wc^3) - (6.93*(1 + 27*upa^3/mat.wc^3)*e^(-6.93*upa/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.σmax0
    end
    return dσmax
end

function calc_kn_ks(mat::MCJoint, ipd::MCJointIpState)
    kn = mat.E*mat.α/ipd.h
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.α/ipd.h
    return kn, ks
end

function find_intersection(mat::MCJoint, ipd::MCJointIpState, F1::Float64, F2::Float64, Δσ::Array{Float64,1})
    @assert(F1*F2<0.0)
    α  = 0.0
    α1 = 0.0
    α2 = 1.0
    F  = F1
    maxits = 40
    for i = 1:maxits
        α  = α1 + F1/(F1-F2)*(α2-α1)
        F  = yield_func(mat, ipd, ipd.σ + α*Δσ)

        if abs(F) < 1e-7 break end

        if F < 0.0
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
    maxits = 50
    Δλ     = 0.0
    f      = 0.0
    upa    = 0.0
    tol    = 1e-2      

    for Δλ in  (0.0, 1e-10)
        for i=1:maxits
            μ      = mat.μ
    		kn, ks = calc_kn_ks(mat, ipd)

		    # quantities at n+1
		    if ipd.ndim == 3
			    if σtr[1]>0
			        σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
			        dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
			        drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
			    else
			        σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
			        dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
			        drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks+1)^2 ]
			    end
			else
			    if σtr[1]>0
			        σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks) ]
			        dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			        drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			    else
			        σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
			        dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			        drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
			    end
			end
			 	
		 	r      = potential_derivs(mat, ipd, σ)
		 	norm_r = norm(r)
		 	upa    = ipd.upa + Δλ*norm_r
		 	σmax   = calc_σmax(mat, ipd, upa)
		    m      = σmax_deriv(mat, ipd, upa)
		    dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

		    if ipd.ndim == 3
		    	f = sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*μ
		    	if (σ[2]==0 && σ[3]==0) 
		        	dfdΔλ = (dσdΔλ[1] - dσmaxdΔλ)*μ		      
		    	else
		        	dfdΔλ = 1/sqrt(σ[2]^2 + σ[3]^2) * (σ[2]*dσdΔλ[2] + σ[3]*dσdΔλ[3]) + (dσdΔλ[1] - dσmaxdΔλ)*μ
		    	end
		    else
				f = abs(σ[2]) + (σ[1]-σmax)*mat.μ
				dfdΔλ = sign(σ[2])*dσdΔλ[2] + (dσdΔλ[1] - dσmaxdΔλ)*μ
		   	end

            Δλ = Δλ - f/dfdΔλ

            abs(f) < tol && break

            if i == maxits || isnan(Δλ)
            	@show i, Δλ
                Δλ = -1.0
                break
            end
        end
        abs(f) < tol && break
    end
    return Δλ
end

function calc_σ_upa(mat::MCJoint, ipd::MCJointIpState, σtr::Array{Float64,1})
    μ = mat.μ
    kn, ks = calc_kn_ks(mat, ipd) 

    if ipd.ndim == 3
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
    ipd.σ = σ
    r = potential_derivs(mat, ipd, ipd.σ)
    ipd.upa += ipd.Δλ*norm(r)
    return ipd.σ, ipd.upa
end

function subincrement(mat::MCJoint, ipd::MCJointIpState, Δwep::Array{Float64,1}, De::Array{Float64,2})
    n = 200
    Δwinc = Δwep/n

    for i = 1:n
        σtr = ipd.σ + De*Δwinc
        ipd.Δλ = calc_Δλ(mat, ipd, σtr)
        if ipd.Δλ < 0 
            if ipd.upa/mat.wc > 0.99
                ipd.upa = mat.wc
                ipd.σ = zeros(ipd.ndim)
                break
            else 
                break 
            end
        end
        if ipd.Δλ > 0; ipd.σ, ipd.upa = calc_σ_upa(mat, ipd, σtr) end
    end
    return ipd.σ, ipd.upa 
end

function mountD(mat::MCJoint, ipd::MCJointIpState)
    kn, ks = calc_kn_ks(mat, ipd)
    σmax = calc_σmax(mat, ipd, ipd.upa)

    if ipd.ndim == 3
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
        v    = yield_deriv(mat, ipd)
        r    = potential_derivs(mat, ipd, ipd.σ)
        y    = -mat.μ # ∂F/∂σmax
        m    = σmax_deriv(mat, ipd, ipd.upa)  # ∂σmax/∂upa

        #Dep  = De - De*r*v'*De/(v'*De*r - y*m*norm(r))

        if ipd.ndim == 3
            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                     -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]
        else
            den = kn*r[1]*v[1] + ks*r[2]*v[2] - y*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  ]
        end

    	return Dep
    end
end

function stress_update(mat::MCJoint, ipd::MCJointIpState, Δw::Array{Float64,1})
    σini = copy(ipd.σ)

    μ = mat.μ
    kn, ks = calc_kn_ks(mat, ipd)
    σmax = calc_σmax(mat, ipd, ipd.upa)  

    if ipd.ndim == 3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn  0.0
               0.0  ks ]
    end

    # σ trial and F trial
    σtr  = ipd.σ + De*Δw
    Ftr  = yield_func(mat, ipd, σtr) 

    # Elastic and EP integration
    if σmax == 0.0 && ipd.w[1] >= 0.0
        if ipd.ndim==3
            r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
            r = r1/norm(r1)
            ipd.Δλ = norm(r1)
        else
            r1 = [ σtr[1]/kn, σtr[2]/ks ]
            r = r1/norm(r1)
            ipd.Δλ = norm(r1)  
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
        Fini = yield_func(mat, ipd, ipd.σ)

        # Elastic increment
        if Ftr > 1e-6 && Fini < 0.0
            α = find_intersection(mat, ipd, Fini, Ftr, De*Δw)
            ipd.w += α*Δw
            ipd.σ += α*De*Δw
        end
  
        σi = copy(ipd.σ)
        upai = copy(ipd.upa)

        # Elastic-plastic increment
        Δwep = (1.0-α)*Δw
        σtr = ipd.σ + De*Δwep
        ipd.Δλ = calc_Δλ(mat, ipd, σtr)
               
        if ipd.Δλ < 0 && ipd.upa/mat.wc > 0.99
        	@show ipd.Δλ, ipd.upa/mat.wc
            ipd.upa = mat.wc
            ipd.σ = zeros(ipd.ndim)
        elseif ipd.Δλ < 0
        	warn("subincrement")
            ipd.σ, ipd.upa = subincrement(mat, ipd, Δwep, De)
        else
            ipd.σ, ipd.upa = calc_σ_upa(mat, ipd, σtr)
                   
           if mat.softcurve == "bilinear"
                if mat.ws > upai  && ipd.upa > mat.ws
                    ipd.σ = σi
                    ipd.upa = upai
                    ipd.σ, ipd.upa = subincrement(mat, ipd, Δwep, De)
                end
            end
        end

        # Return to surface:
        F  = yield_func(mat, ipd, ipd.σ)

     
        #@show ipd.Δλ
        #if false 
        #if F > 1e-1
        #    ipd.σ = σi
        #    ipd.upa = upai
        #    ipd.σ, ipd.upa = subincrement(mat, ipd, Δwep, De)
        #end
  

        #  Plastic update of w and Δσ
        ipd.w += Δwep
        Δσ = ipd.σ - σini
    end
    F  = yield_func(mat, ipd, ipd.σ)

    if F > 1e-1
    	@show F
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
          :upa => ipd.upa
          )
    else
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :upa => ipd.upa
          )
    end
end
