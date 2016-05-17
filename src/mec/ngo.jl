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


export NgoScordelis

type NgoScordelisIpData<:IpData
    ndim::Int
    sig ::Array{Float64,1}
    eps ::Array{Float64,1}
    function NgoScordelisIpData(ndim=3)
        this = new(ndim)
        this.sig = zeros(3)
        this.eps = zeros(3)
        return this
    end
end

type NgoScordelis<:AbsJoint
    ks::Float64
    kn::Float64
    new_ipdata::DataType

    function NgoScordelis(prms::Dict{Symbol,Float64})
        return  NgoScordelis(;prms...)
    end

    function NgoScordelis(;ks=NaN, kn=NaN, ft=0.0)
        @assert ks>=0
        @assert kn>=0
        @assert ft>=0

        this = new(ks, kn)
        this.new_ipdata = NgoScordelisIpData
        return this
    end
end

function set_state(ipd::NgoScordelisIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.sig[:] = sig
    else
        if length(sig)!=0; error("NgoScordelis: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.eps[:] = eps
    else
        if length(eps)!=0; error("NgoScordelis: Wrong size for strain array: $eps") end
    end
end

function Tau(mat::NgoScordelis, s::Float64)
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

function Tau_deriv(mat::NgoScordelis, ipd::CEBJoint1DIpData, s::Float64)
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

function Sig(mat::NgoScordelis, w::Float64)
    ww = abs(w)
end

function Sig_deriv(mat::NgoScordelis, w::Float64)
    ww = abs(w)
end

function mountD(mat::NgoScordelis, ipd::NgoScordelisIpData)


    s  = abs(ipd.ε[1])
    kh = deriv(mat, ipd, s)
    ks = ipd.unload? mat.ks : mat.ks*kh/(mat.ks + kh)

    kn = mat.kn
    if ipd.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end



    ks = mat.ks
    kn = mat.kn
    if ipd.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   kn ]
    end
end

function stress_update(mat::NgoScordelis, ipd::NgoScordelisIpData, Δu)
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.eps[1:ipd.ndim] += Δu
    ipd.sig[1:ipd.ndim] += Δσ
    return Δσ
end

function getvals(mat::NgoScordelis, ipd::NgoScordelisIpData)
    if ipd.ndim == 2
        return Dict(
          :s11  => ipd.sig[1] ,
          :s12  => ipd.sig[2] )
    else
        return Dict(
          :s11  => ipd.sig[1] ,
          :s12  => ipd.sig[2] ,
          :s13  => ipd.sig[3] )
    end
end

function node_and_elem_vals(mat::NgoScordelis, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    return node_vals, elem_vals
end

