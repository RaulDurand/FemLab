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


export NLJoint

# Integration point data
type NLJointIpData<:IpData
    ndim::Int
    σ ::Array{Float64,1} # stress vector (τx', τy', σn)
    u ::Array{Float64,1} # relative displacement (sx', sy', w)
    h ::Float64          # surrouding element size h=(V1+V2)/(2A)
    function NLJointIpData(ndim=3)
        this = new(ndim)
        this.σ = zeros(3)
        this.u = zeros(3)
        this.h = 0.0
        return this
    end
end

# Material
type NLJoint<:AbsJoint
    fc::Float64  # surrounding material compressive strength
    ft::Float64  # surrounding material tensile strength
    wc::Float64  # critical crack opennig (e.g. 160 μm for concrete)
    new_ipdata::DataType

    function NLJoint(prms::Dict{Symbol,Float64})
        return  NLJoint(;prms...)
    end

    function NLJoint(;fc=NaN, ft=NaN, wc=NaN)
        @assert fc >= 0.0
        @assert wc >  0.0
        if isnan(ft)
            ft = 0.1*fc
        end
        @assert ft >= 0.0

        this = new(fc, ft, wc)
        this.new_ipdata = NLJointIpData
        return this
    end
end

const CWT = 0.005
const CWB = 0.3
const CFB = 0.1
const CFM = 0.5

function set_state(ipd::NLJointIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.σ[:] = sig
    else
        if length(sig)!=0; error("NLJoint: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.u[:] = eps
    else
        if length(eps)!=0; error("NLJoint: Wrong size for strain array: $eps") end
    end
end

function Tau(mat::NLJoint, ipd::NLJointIpData, s::Float64, w::Float64)
    ss = abs(s)
    wt = mat.wc*CWT.*ipd.h
    fm = CFM*mat.ft

    #if w>wt return 0.0 end

    #if ss<wt
        ks = fm/wt
        ks = 3e+6
        return ks*s
    #end

    wc = mat.wc
    fb = CFB*fm
    wb = wt + CWB*(wc-wt)
    if ss<wb
        return (fm + (fb-fm)/(wb-wt)*(ss-wt))*sign(s)
    elseif ss<wc
        return (fb - fb/(wc-wb)*(ss-wb))*sign(s)
    else
        return 0.0
    end
end

function Tau_deriv(mat::NLJoint, ipd::NLJointIpData, s::Float64, w::Float64)
    #@show s
    ss = abs(s)
    wt = mat.wc*CWT.*ipd.h
    fm = CFM*mat.ft

    #if w>wt return 0.0 end

    #if ss<wt
        ks = fm/wt
        ks = 3e+6
        return ks
    #end

    wc = mat.wc
    fb = CFB*fm
    wb = wt + CWB*(wc-wt)
    if ss<wb
        return (fb-fm)/(wb-wt)
    elseif ss<wc
        return -fb/(wc-wb)
    else
        return 0.0
    end
end

function Sig(mat::NLJoint, ipd::NLJointIpData, w::Float64)
    #wt = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)
    wt = mat.wc*CWT.*ipd.h

    if w<0.0
        kn = 10.0*mat.ft/wt
        return kn*w
    elseif w<wt
        kn = mat.ft/wt
        return kn*w
    end

    wc = mat.wc
    fm = mat.ft
    fb = CFB*fm
    #wb = 0.5*(wt+wc)
    wb = wt + CWB*(wc-wt)
    if w<wb
        return fm + (fb-fm)/(wb-wt)*(w-wt)
    elseif w<wc
        return fb - fb/(wc-wb)*(w-wb)  
    else
        return 0.0
    end
end

function Sig_deriv(mat::NLJoint, ipd::NLJointIpData, w::Float64)

    wt = mat.wc*CWT.*ipd.h

    if w<0.0
        kn = 10.0*mat.ft/wt
        return kn
    elseif w<wt
        kn = mat.ft/wt
        return kn
    end

    wc = mat.wc
    fm = mat.ft
    fb = CFB*fm
    #wb = 0.5*(wt+wc)
    wb = wt + CWB*(wc-wt)
    if w<wb
        return (fb-fm)/(wb-wt)
    elseif w<wc
        return -fb/(wc-wb)
    else
        #return 1e-3*fm/wt
        return 0.0
    end
end


function mountD(mat::NLJoint, ipd::NLJointIpData)
    s1 = ipd.u[1]
    s2 = ipd.u[2]
    w  = ipd.u[3]

    ks1 = Tau_deriv(mat, ipd, s1, w)
    ks2 = Tau_deriv(mat, ipd, s2, w)
    kn  = Sig_deriv(mat, ipd, w)
    #@show kn
    #@show ks1
    #exit()
    #@show wt
    #@show w, s1, s2
    #@show kn, ks1, ks2

    if ipd.ndim==2
        return [ ks1  0.0 
                 0.0  kn ]
    else
        return [ ks1  0.0  0.0
                 0.0  ks2  0.0
                 0.0  0.0  kn ]
    end
end

function stress_update(mat::NLJoint, ipd::NLJointIpData, Δu)
    s1 = ipd.u[1]
    s2 = ipd.u[2]
    w  = ipd.u[3]

    # relative displacements
    Δs1 = Δu[1]
    Δs2 = Δu[2]
    Δw  = Δu[3]
    #@show Δw

    # increment of stress
    Δσ = zeros(3)

    τ1 = Tau(mat, ipd, s1+Δs1, w)
    τ2 = Tau(mat, ipd, s2+Δs2, w)
    σn = Sig(mat, ipd, w+Δw)
    #@show w+Δw
    #@show σn

    Δσ[1] = τ1 - ipd.σ[1]
    Δσ[2] = τ2 - ipd.σ[2]
    Δσ[3] = σn - ipd.σ[3]

    # update of relative displacements and stresses
    ipd.u[1:ipd.ndim] += Δu
    ipd.σ[1:ipd.ndim] += Δσ

    return Δσ
end

function getvals(mat::NLJoint, ipd::NLJointIpData)
    #wt  = 2.0*ipd.h*mat.ft/(mat.p0*mat.E)

    if ipd.ndim == 2
        return Dict(
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] )
    else
        w      = ipd.u[3]
        #ncrack = w>wt ? 1.0: 0.0
        return Dict(
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :sn  => ipd.σ[3] ,
          :w    => ipd.u[3] ,
          #:wt   => wt ,
          #:ncrack => ncrack
          )
    end
end

function node_and_elem_vals(mat::NLJoint, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    return node_vals, elem_vals
end

