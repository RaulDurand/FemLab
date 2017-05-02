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


export Joint1D

type Joint1DIpData<:IpData
    ndim::Int
    sig ::Array{Float64,1}
    eps ::Array{Float64,1}
    function Joint1DIpData(ndim=3)
        this = new(ndim)
        this.sig = zeros(3)
        this.eps = zeros(3)
        return this
    end
end

type Joint1D<:AbsJoint1D
    ks::Float64
    kn::Float64
    h ::Float64    # section perimeter

    function Joint1D(prms::Dict{Symbol,Float64})
        return  Joint1D(;prms...)
    end

    function Joint1D(;ks=NaN, kn=NaN, h=NaN, A=NaN, dm=NaN)
        # A : section area
        # dm: section diameter
        # h : section perimeter
        @assert ks>=0
        @assert kn>=0
        @assert (h>0 || A>0 || dm>0)

        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^0.5
            else
                h = pi*dm
            end
        end
        @assert h>0

        this = new(ks, kn, h)
        return this
    end
end

# Create a new instance of Ip data
new_ipdata(mat::Joint1D, ndim::Int) = Joint1DIpData(ndim)

function set_state(ipd::Joint1DIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.sig[:] = sig
    else
        if length(sig)!=0; error("MecElasticSolid: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.eps[:] = eps
    else
        if length(eps)!=0; error("MecElasticSolid: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::Joint1D, ipd::Joint1DIpData)
    ks = mat.ks
    kn = mat.kn
    if ipd.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function stress_update(mat::Joint1D, ipd::Joint1DIpData, deps)
    D = calcD(mat, ipd)
    dsig = D*deps

    ipd.eps[1:ipd.ndim] += deps
    ipd.sig[1:ipd.ndim] += dsig
    return dsig
end

function getvals(mat::Joint1D, ipd::Joint1DIpData)
    return Dict(
      :ur   => ipd.eps[1] ,
      :tau  => ipd.sig[1] )
end
