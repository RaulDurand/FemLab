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


export Hordijk

# Integration point data
type HordijkIpData<:IpData
    ndim::Int
    sig ::Array{Float64,1} # stress vector (τx', τy', σn)
    eps ::Array{Float64,1} # relative displacement (sx', sy', w)
    h   ::Float64          # surrouding element size h=(V1+V2)/(2A)
    function HordijkIpData(ndim=3)
        this = new(ndim)
        this.sig = zeros(3)
        this.eps = zeros(3)
        this.h = 0.0
        return this
    end
end

# Material
type Hordijk<:AbsJoint
    ks::Float64  # joint shear stiffness
    E ::Float64  # surrounding material Young's modulus
    ft::Float64  # surrounding material tensile strength
    wc::Float64  # critical crack opennig (e.g. 160 μm for concrete)
    p0::Float64  # penalty value
    c ::Float64  # constant to allow curve smoothing
    new_ipdata::DataType

    function Hordijk(prms::Dict{Symbol,Float64})
        return  Hordijk(;prms...)
    end

    function Hordijk(;ks=NaN, E=NaN, fc=NaN, ft=NaN, wc=NaN, p0=100.0, c=1.0)
        @assert ks>=0

        if fc>0.0
            #ft = funcao... TODO
        end
        @assert E  >= 0.0
        @assert ft >= 0.0
        @assert wc >  0.0
        @assert p0 >= 0.0
        @assert c >= 1.0 && c<=2.0

        this = new(ks, E, ft, wc, p0, c)
        this.new_ipdata = HordijkIpData
        return this
    end
end

function set_state(ipd::HordijkIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.sig[:] = sig
    else
        if length(sig)!=0; error("Hordijk: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.eps[:] = eps
    else
        if length(eps)!=0; error("Hordijk: Wrong size for strain array: $eps") end
    end
end

function z(ρ, c)
    ρ0   = ρ^c
    return (1.0+27.0*ρ0^3.0)*exp(-7.0*ρ0) - 28.0*ρ0*exp(-7.0)
    #return (1.0+27.0*ρ^3.0)*exp(-7.0*ρ) - 28.0*ρ*exp(-7.0)
end

function Tau(mat::Hordijk, ipd::HordijkIpData, s::Float64, w::Float64)
    wt = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)
    st = wt

    if w<wt
        return 2*s/st*mat.ft   # z==1.0
    elseif w<mat.wc
        ρ = (w-wt)/(mat.wc-wt) # varies from 0 to 1
        c = mat.c
        return 2.0*s/st*mat.ft*z(ρ,c)
    else 
        return 0.0
    end
end

function Tau_deriv(mat::Hordijk, ipd::HordijkIpData, s::Float64, w::Float64)
    wt = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)
    st = wt
    if w<wt
        return 0.0, 2/st*mat.ft  # dτ/dw dτ/ds
    elseif w<mat.wc
        ρ    = (w-wt)/(mat.wc-wt) # varies from 0 to 1
        dρdw = 1.0/(mat.wc - wt)
        c    = mat.c
        ρ0   = ρ^c
        #dzdρ = 81.0*ρ^2.0*exp(-7.0*ρ) - 7.0*exp(-7.0*ρ) - 189.0*(ρ^3.0)*exp(-7.0*ρ) - 28.0*exp(-7.0)
        dzdρ0 = 81.0*ρ0^2.0*exp(-7.0*ρ0) - 7.0*exp(-7.0*ρ0) - 189.0*(ρ0^3.0)*exp(-7.0*ρ0) - 28.0*exp(-7.0)
        dρ0dρ = c*ρ^(c-1)
        dzdρ  = dzdρ0*dρ0dρ
        dzdw  = dzdρ*dρdw
        return 2*abs(s)/st*mat.ft*dzdw, 2/st*mat.ft*z(ρ,c)  # dτ/dw dτ/ds
    else
        return 0.0, 0.0
    end
end

function Sig(mat::Hordijk, ipd::HordijkIpData, w::Float64)
    wt = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)
    if w<0.0
        return 2*w/wt*mat.ft
    elseif w<wt
        kn0 = mat.ft/wt
        return kn0*w
    elseif w<mat.wc
        ρ = (w-wt)/(mat.wc-wt) # varies from 0 to 1
        c = mat.c
        return z(ρ,c)*mat.ft
    else
        return 0.0
    end
end

function Sig_deriv(mat::Hordijk, ipd::HordijkIpData, w::Float64)
    wt = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)
    kn0 = mat.ft/wt
    if w<0.0
        return 2*kn0
    elseif w<wt
        return kn0
    elseif w<mat.wc
        # dσ/dw = dz/dw*ft
        # dz/dw = dz/dρ*dρ/dw
        dρdw = 1.0/(mat.wc - wt)
        ρ    = (w-wt)/(mat.wc-wt) # varies from 0 to 1
        c    = mat.c
        ρ0   = ρ^c
        dzdρ0 = 81.0*ρ0^2.0*exp(-7.0*ρ0) - 7.0*exp(-7.0*ρ0) - 189.0*(ρ0^3.0)*exp(-7.0*ρ0) - 28.0*exp(-7.0)
        #dzdρ = 81.0*ρ^2.0*exp(-7.0*ρ) - 7.0*exp(-7.0*ρ) - 189.0*(ρ^3.0)*exp(-7.0*ρ) - 28.0*exp(-7.0)
        dρ0dρ = c*ρ^(c-1)
        dzdρ  = dzdρ0*dρ0dρ
        dzdw  = dzdρ*dρdw
        dσdw  = dzdw*mat.ft

        # Tamming the derivative to avoid convergence problems for high values
        #max_dσdw = kn0/10000.0
        #if abs(dσdw) > max_dσdw
            #dσdw = abs(dσdw)/dσdw * max_dσdw
        #end

        return dσdw
    else
        return 0.0
    end
end

function mountD(mat::Hordijk, ipd::HordijkIpData)

    wt = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)
    ks = mat.ks
    w  = ipd.eps[3]

    kn0 = mat.ft/wt
    kn  = Sig_deriv(mat, ipd, w)

    if w>wt
        kn = kn0*kn/(kn0+kn)
    end

    #if ipd.ndim==2
        #return [  ks  0.0 
                 #0.0   kn ]
    #else
        #return  [  ks  0.0  0.0
                  #0.0   ks  0.0
                  #0.0  0.0   kn ]
    #end

    if ipd.ndim==2
        #  D = [ dτ1/ds1  dτ1/dw ]
        #      [     0.0   dσ/dw ]

        s1 = ipd.eps[1]
        ks1, dτ1dw = Tau_deriv(mat, ipd, s1, w)
        return [   ks1   dτ1dw 
                   0.0      kn ]
    else
        s1 = ipd.eps[1]
        s2 = ipd.eps[2]
        dτ1dw, ks1 = Tau_deriv(mat, ipd, s1, w)
        dτ2dw, ks2 = Tau_deriv(mat, ipd, s2, w)
        #@show ks1
        #@show ks2
        #@show dτ1dw
        #@show dτ2dw

        #      [ dτ1/ds1      0.0  dτ1/dw ]
        #  D = [     0.0  dτ2/ds2  dτ2/dw ]
        #      [     0.0      0.0   dσ/dw ]

        return [   ks1   0.0  dτ1dw 
                   0.0   ks2  dτ2dw 
                   0.0   0.0     kn ]

        #return [   ks1   0.0    0.0
                   #0.0   ks2    0.0
                   #0.0   0.0     kn ]
        #ks = mat.ks
        #return [   ks    0.0    0.0
                   #0.0   ks     0.0
                   #0.0   0.0     kn ]
    end
end

function stress_update(mat::Hordijk, ipd::HordijkIpData, Δu)
    #@show Δu
    wt  = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)
    @assert wt < mat.wc

    kn0 = mat.ft/wt
    ks  = mat.ks

    # relative displacements
    Δw = Δu[3]

    # increment of relative displacements
    Δσ = zeros(3)

    # relative displacements at end of increment
    wf = ipd.eps[3] + Δw
    s1f = ipd.eps[1] + Δu[1]
    s2f = ipd.eps[2] + Δu[2]

    # normal stress increment
    if wf<0.0
        Δσ[3] = 2*kn0*Δw
    elseif wf<wt
        Δσ[3] = kn0*Δw
    else
        σn = Sig(mat, ipd, wf)
        Δσ[3] = σn - ipd.sig[3]
    end

    # shear stress increments
    Δσ[1] = ks*Δu[1]
    Δσ[2] = ks*Δu[2]
    if wf<wt
        st = wt
        ks1 = ks2 = 2/st*mat.ft
        Δσ[1] = ks1*Δu[1]
        Δσ[2] = ks2*Δu[2]
    else
        τ1 = Tau(mat, ipd, s1f, wf) 
        Δσ[1] = τ1 - ipd.sig[1]
        τ2 = Tau(mat, ipd, s2f, wf) 
        Δσ[2] = τ2 - ipd.sig[2]
    end

    # update of relative displacements and stresses
    ipd.eps[1:ipd.ndim] += Δu
    #@show ipd.eps

    ipd.sig[1:ipd.ndim] += Δσ

    return Δσ
end

function getvals(mat::Hordijk, ipd::HordijkIpData)
    wt  = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)

    if ipd.ndim == 2
        return Dict(
          :s1  => ipd.sig[1] ,
          :s2  => ipd.sig[2] )
    else
        w      = ipd.eps[3]
        ncrack = w>wt ? 1.0: 0.0

        s1 = ipd.eps[1]
        dτ1dw, ks1 = Tau_deriv(mat, ipd, s1, w)
        return Dict(
          :s1  => ipd.eps[1] ,
          :s2  => ipd.eps[2] ,
          :w    => ipd.eps[3] ,
          :wt   => wt ,
          :tau1 => ipd.sig[1] ,
          :tau2 => ipd.sig[2] ,
          :s_n  => ipd.sig[3] ,
          :ks1 => ks1,
          :dt1dw => dτ1dw,
          :ncrack => ncrack
          )
    end
end

function node_and_elem_vals(mat::Hordijk, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    return node_vals, elem_vals
end

