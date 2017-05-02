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


export Joint

type JointIpData<:IpData
    ndim::Int
    σ   ::Array{Float64,1}
    w   ::Array{Float64,1}
    h   ::Float64
    function JointIpData(ndim=3)
        this = new(ndim)
        this.σ   = zeros(3)
        this.w = zeros(3)
        this.h = 0.0
        return this
    end
end

type Joint<:AbsJoint
    E::Float64
    ν::Float64
    α::Float64

    function Joint(prms::Dict{Symbol,Float64})
        return  Joint(;prms...)
    end

    function Joint(;E=NaN, nu=NaN, alpha=1.0)
        @assert E>=0
        @assert nu>=0

        this = new(E, nu, alpha)
        return this
    end
end

# Create a new instance of Ip data
new_ipdata(mat::Joint, ndim::Int) = JointIpData(ndim)

function set_state(ipd::JointIpData, sig=zeros(0), w=zeros(0))
    if length(sig)==3
        ipd.σ[:] = sig
    else
        if length(sig)!=0; error("Joint: Wrong size for stress array: $sig") end
    end
    if length(w)==3
        ipd.w[:] = w
    else
        if length(w)!=0; error("Joint: Wrong size for strain array: $w") end
    end
end

function mountD(mat::Joint, ipd::JointIpData)
    G  = mat.E/(1.0+mat.ν)/2.0
    kn = mat.E*mat.α/ipd.h
    ks =     G*mat.α/ipd.h
    if ipd.ndim==2
        return [  kn  0.0 
                 0.0   ks ]
    else
        return  [  kn  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   ks ]
    end
end

function stress_update(mat::Joint, ipd::JointIpData, Δu)
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.w[1:ipd.ndim] += Δu
    ipd.σ[1:ipd.ndim] += Δσ
    return Δσ
end

function getvals(mat::Joint, ipd::JointIpData)
    if ipd.ndim == 2
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] )
    else
        return Dict(
          :w1  => ipd.w[1] ,
          :w2  => ipd.w[2] ,
          :w3  => ipd.w[3] ,
          :s1  => ipd.σ[1] ,
          :s2  => ipd.σ[2] ,
          :s3  => ipd.σ[3] )
    end
end
