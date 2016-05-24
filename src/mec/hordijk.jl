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
    new_ipdata::DataType

    function Hordijk(prms::Dict{Symbol,Float64})
        return  Hordijk(;prms...)
    end

    function Hordijk(;ks=NaN, E=NaN, fc=NaN, ft=NaN, wc=NaN)
        @assert ks>=0

        if fc>0.0
            #ft = funcao... TODO
        end
        @assert E>=0
        @assert ft>=0
        @assert wc>0

        this = new(ks, E, ft, wc)
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

function Tau(mat::Hordijk, s::Float64)
    ss = abs(s)
    #if ss<mat.s1
        #return mat.ks*s
    #elseif ss<mat.s2
        #return mat.ks*mat.s1*sign(s)
    #elseif ss<mat.s3
        #τmax = mat.ks*mat.s1
        #return (mat.τres + (mat.τres-τmax)/(mat.s3-mat.s2)*(ss-mat.s3))*sign(s)
        #
    #else
        #return mat.τres*sign(s)
    #end
end

function Tau_deriv(mat::Hordijk, ipd::CEBJoint1DIpData, s::Float64)
    ss = abs(s)
    
    #if ss<mat.s1
        #return mat.ks
    #elseif ss<mat.s2
        #return 0.0
    #elseif ss<mat.s3
        #τmax = mat.ks*mat.s1
        #return (mat.τres-τmax)/(mat.s3-mat.s2)
    #else
        #return 0.0
    #end
end

function Sig(mat::Hordijk, ipd::HordijkIpData, w::Float64)
    wt = 2.0*ipd.h*mat.ft/(100.0*mat.E)
    if w<wt
        kn0 = mat.ft/wt
        return kn0*w
    elseif w<mat.wc
        ρ = (w-wt)/(mat.wc-wt) # varies from 0 to 1
        z = (1.0+27.0*ρ^3.0)*exp(-7.0*ρ) - 28.0*ρ*exp(-7.0)
        return z*mat.ft
    else
        return 0.0
    end
end

function Sig_deriv(mat::Hordijk, ipd::HordijkIpData, w::Float64)
    wt = 2.0*ipd.h*mat.ft/(100.0*mat.E)
    if w<wt
        kn0 = mat.ft/wt
        return kn0
    elseif w<mat.wc
        # dσ/dw = dz/dw*ft
        # dz/dw = dz/dρ*dρ/dw
        dρdw = 1.0/(mat.wc - wt)
        ρ    = (w-wt)/(mat.wc-wt) # varies from 0 to 1
        dzdρ = 81.0*ρ^2.0*exp(-7.0*ρ) - 7*exp(-7.0*ρ) - 189.0*(ρ^3)*exp(-7.0*ρ) - 28.0*exp(-7.0)
        dzdw = dzdρ*dρdw
        return dzdw*mat.ft
    else
        return 0.0
    end
end

function mountD(mat::Hordijk, ipd::HordijkIpData)

    wt = 2.0*ipd.h*mat.ft/(100.0*mat.E)
    ks = mat.ks
    w  = ipd.eps[3]

    kn0 = mat.ft/wt
    kn  = Sig_deriv(mat, ipd, w)

    if w>wt
        kn = kn0*kn/(kn0+kn)
    end

    if ipd.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   kn ]
    end
end

function stress_update(mat::Hordijk, ipd::HordijkIpData, Δu)
    wt  = 2.0*ipd.h*mat.ft/(100.0*mat.E)
    kn0 = mat.ft/wt
    ks  = mat.ks

    # relative displacements
    s1 = Δu[1]
    s2 = Δu[2]
    Δw = Δu[3]

    # increment of relative displacements
    Δσ = zeros(3)

    # relative displacements at end of increment
    wf = ipd.eps[3] + Δw

    # normal stress increment
    if wf<wt
        Δσ[3] = kn0*Δw
    else
        σn = Sig(mat, ipd, wf)
        Δσ[3] = σn - ipd.sig[3]
    end

    # shear stress increments
    Δσ[1] = ks*Δu[1]
    Δσ[2] = ks*Δu[2]

    # update of relative displacements and stresses
    ipd.eps[1:ipd.ndim] += Δu
    ipd.sig[1:ipd.ndim] += Δσ
    return Δσ
end

function getvals(mat::Hordijk, ipd::HordijkIpData)
    if ipd.ndim == 2
        return Dict(
          :s1  => ipd.sig[1] ,
          :s2  => ipd.sig[2] )
    else
        return Dict(
          :s1  => ipd.sig[1] ,
          :s2  => ipd.sig[2] ,
          :s_n  => ipd.sig[3] ,
          :w    => ipd.eps[3] ,
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

