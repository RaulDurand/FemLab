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
        b = 1e-5
    end
    return a, b
end

function yield_func(mat::MCJoint, ipd::MCJointIpData, σ::Array{Float64,1}, upa)
    a, b = calc_σmax(mat, upa)
    σmax = a - b*upa
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
        r = [0.0, 2*σ[2], 2*σ[3]]
        return r/norm(r)
    else
        τmax = σmax*mat.μ
        r = [2*σ[1]/σmax^2, 2*σ[2]/τmax^2, 2*σ[3]/τmax^2]
        return r/norm(r)
        #return [2*σ[1]/σmax^2, 2*σ[2]/τmax^2, 2*σ[3]/τmax^2]
    end
end

function normβ(V::Array{Float64,1}, β::Float64)
    return (V[1]^2 + β*V[2]^2 + β*V[3]^2)^0.5
end

function mountD(mat::MCJoint, ipd::MCJointIpData)
    kn = mat.E*mat.α/ipd.h
    #if ipd.w[1]<0.0
        #kn = mat.E*mat.α/ipd.h*10.0
    #end

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
        #@show Dep
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
    #@show kn*mat.ws
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
    else
        local a, b, r, nr, q
        
        eps = 1e-7
        err = eps + 1
        x0  = 1.0
        nmaxits = 40
        f = 0
        i=0
        x0 = ipd.upa

        function func(x)
            # at n
            r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)  # r_n
            nr= normβ(r, mat.β)
            Δupa = x*nr
            q = De*r
            s   = σtr - x*q

            # at n+1
            r = potential_derivs(mat, ipd, s, ipd.upa+Δupa)  # r_n+1 
            nr= normβ(r, mat.β)
            Δupa = x*nr
            q = De*r
            s   = σtr - x*q

            # at n+1
            r = potential_derivs(mat, ipd, s, ipd.upa+Δupa)  # r_n+1 
            nr= normβ(r, mat.β)
            Δupa = x*nr
            q = De*r
            s   = σtr - x*q

            a, b = calc_σmax(mat, ipd.upa+x*nr)
            root = ( (σtr[2] - x*q[2])^2 + (σtr[3] - x*q[3])^2 )^0.5 
            f = root + (σtr[1] - x*q[1])*mat.μ - (a - b*(ipd.upa + x*nr))*mat.μ
            d = ( (x*q[2] - σtr[2])*q[2] + (x*q[3] - σtr[3])*q[3] )/root + (b*nr -q[1])*mat.μ
            return f, d, r, nr, q, a, b
        end


        change_deriv_sgn = false
        xx = 0.0
        for i=1:40

            #if false
            # at n
            r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)  # r_n
            nr= normβ(r, mat.β)
            Δupa = x0*nr
            q = De*r
            s   = σtr - x0*q

            # at n+1
            r = potential_derivs(mat, ipd, s, ipd.upa+Δupa)  # r_n+1 
            nr= normβ(r, mat.β)
            Δupa = x0*nr
            q = De*r
            s   = σtr - x0*q

            # at n+1
            r = potential_derivs(mat, ipd, s, ipd.upa+Δupa)  # r_n+1 
            nr= normβ(r, mat.β)
            Δupa = x0*nr
            q = De*r
            s   = σtr - x0*q


            a, b = calc_σmax(mat, ipd.upa+x0*nr)
            root = ( (σtr[2] - x0*q[2])^2 + (σtr[3] - x0*q[3])^2 )^0.5 
            f = root + (σtr[1] - x0*q[1])*mat.μ - (a - b*(ipd.upa + x0*nr))*mat.μ
            d = ( (x0*q[2] - σtr[2])*q[2] + (x0*q[3] - σtr[3])*q[3] )/root + (b*nr -q[1])*mat.μ
            #end


            #f, d, r, nr, q, a, b = func(x0)

            if change_deriv_sgn
                if d>0 && a>0
                    d = -d*10
                end
                #else
                #d = -abs(d)*10.0
            end

            x = x0 - f/d

            err = abs(x- x0)

            x0 = x
            if err<eps 
                break
                if x0>0.0 break end
                #@show i,x0
                xx = x
                x0 = 0.0
                #@show i,x0
                change_deriv_sgn = true
            end
        end

        if x0<0.0
            x0 = 1e-10
            f, d, r, nr, q, a, b = func(x0)
        end

        ipd.Δλ = x0

        #if i==40 && xx!=0.0
            #ipd.Δλ = xx
        #end


        #@assert(ipd.Δλ>0)
        #if ipd.Δλ < 0.0
            #pcolor(:red, "Warning: ipd.Δλ<0 ($(ipd.Δλ)) \n")
        #end


        if false
        if ipd.Δλ<0 || i== 40
            @show i
            @show yield_func(mat, ipd, σini, ipd.upa) 
            @show r
            @show ipd.Δλ
            @show σini
            @show σtr
            @show Δw
            @show a,b
            @show ipd.upa + x0*nr
            @show err
            @show yield_func(mat, ipd, σtr - x0*q, ipd.upa + x0*nr) 
            X = linspace(-0.001,0.001,1000)
            Y = []
            F = []


            for x in X
                #s    = σtr - x*q
                r = potential_derivs(mat, ipd, ipd.σ, ipd.upa)  # r_n+1 -> r_n
                nr= normβ(r, mat.β)
                Δupa = x*nr
                q = De*r
                s   = σtr - x*q

                r = potential_derivs(mat, ipd, s, ipd.upa+Δupa)  # r_n+1 -> r_n
                nr= normβ(r, mat.β)
                Δupa = x*nr
                q = De*r
                s   = σtr - x*q

                r = potential_derivs(mat, ipd, s, ipd.upa+Δupa)  # r_n+1 -> r_n
                nr= normβ(r, mat.β)
                Δupa = x*nr
                q = De*r
                s   = σtr - x*q

                r = potential_derivs(mat, ipd, s, ipd.upa+Δupa)  # r_n+1 -> r_n
                nr= normβ(r, mat.β)
                Δupa = x*nr
                q = De*r
                s   = σtr - x*q

                y = yield_func(mat, ipd, s, ipd.upa + Δupa) 
                root = ( (σtr[2] - x*q[2])^2 + (σtr[3] - x*q[3])^2 )^0.5 
                a, b = calc_σmax(mat, ipd.upa+Δupa)
                #@show σini
                #@show a,b
                #@show r
                #@show s
                f = root + (σtr[1] - x*q[1])*mat.μ - (a - b*(ipd.upa + x*nr))*mat.μ
                push!(Y, y)
                push!(F, f)
            end

            PyPlot.plot([-0.001,0.001],[0,0], "-")
            PyPlot.plot(X,Y, "-o")
            PyPlot.plot(X,F, "-^")
            show()
            exit()
        end
        end


        ipd.σ   = σtr - ipd.Δλ*q
        Δupa    = ipd.Δλ*nr
        ipd.upa += Δupa
        #if ipd.upa<0.0
            #ipd.upa = 0.0
            #pcolor(:red, "Warning: ipd.upa<0 ($(ipd.upa))  σmax=$(a-b*ipd.upa)\n")
        #end
        #f = yield_func(mat, ipd, ipd.σ, ipd.upa) 
        #if abs(f)>10
            #pcolor(:red, "Warning: yield function f=$f\n")
        #end
    end
    #@show ipd.upa

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

