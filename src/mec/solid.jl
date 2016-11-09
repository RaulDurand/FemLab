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

export ElasticSolid
export set_state

type ElasticSolidIpData<:IpData
    ndim::Int
    σ::Array{Float64,1}
    ε::Array{Float64,1}
    function ElasticSolidIpData(ndim=3) 
        this = new(ndim)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this
    end
end

type ElasticSolid<:Mechanical
    E ::Float64
    nu::Float64
    De::Tensor4
    new_ipdata::DataType

    function ElasticSolid(prms::Dict{Symbol,Float64})
        return  ElasticSolid(;prms...)
    end

    function ElasticSolid(;E=1.0, nu=0.0)
        if E<=0.0      ; error("Invalid value for E: $E") end
        if !(0<=nu<0.5); error("Invalid value for nu: $nu") end
        this    = new(E, nu)
        this.De = zeros(6,6)
        setDe(E, nu, this.De)
        this.new_ipdata = ElasticSolidIpData

        return this
    end
end

function set_state(ipd::ElasticSolidIpData; sig=zeros(0), eps=zeros(0))
    sq2 = 2.0^0.5
    mdl = [1, 1, 1, sq2, sq2, sq2]
    if length(sig)==6
        ipd.σ[:] = sig.*mdl
    else
        if length(sig)!=0; error("ElasticSolid: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*mdl
    else
        if length(eps)!=0; error("ElasticSolid: Wrong size for strain array: $eps") end
    end
end

function setDe(E::Float64, nu::Float64, De::Array{Float64,2})
    c = E/((1.0+nu)*(1.0-2.0*nu))
    De[:] = 0.0
    De[1,1] = De[2,2] = De[3,3] = c*(1.-nu)
    De[1,2] = De[1,3] = De[2,1] = De[2,3] = De[3,1] = De[3,2] = c*nu
    De[4,4] = De[5,5] = De[6,6] = c*(1.-2.*nu)
end

function calcD(mat::ElasticSolid, ipd::ElasticSolidIpData)
    return mat.De
end

function stress_update(mat::ElasticSolid, ipd::ElasticSolidIpData, dε::Array{Float64,1})
    dσ = mat.De*dε
    ipd.ε += dε
    ipd.σ += dσ
    return dσ
end

function getvals(mat::ElasticSolid, ipd::ElasticSolidIpData)
    σ  = ipd.σ
    ε  = ipd.ε
    ndim = ipd.ndim
    sr2  = 2.0^0.5

    if ndim==2;
        return Dict(
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :sxy => σ[4]/sr2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :exy => ε[4]/sr2,
          :s_m => sum(σ[1:3])/3.0 )
      else
        return Dict(
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :sxy => σ[4]/sr2,
          :syz => σ[5]/sr2,
          :sxz => σ[6]/sr2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :exy => ε[4]/sr2,
          :eyz => ε[5]/sr2,
          :exz => ε[6]/sr2,
          :s_m => sum(σ[1:3])/3.0 )
      end
end
