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
    new_ipdata::DataType

    function ElasticSolid(;E=1.0, nu=0.0)
        if E<=0.0      ; error("Invalid value for E: $E") end
        if !(0<=nu<0.5); error("Invalid value for nu: $nu") end
        this = new(E, nu)
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

function stress_update(mat::ElasticSolid, ipd::ElasticSolidIpData, deps::Array{Float64,1})
    D    = mount_D(mat, ipd)
    dsig = D*deps
    ipd.ε += deps
    ipd.σ += dsig
    dsig
end

function mount_De(E::Float64, nu::Float64)
    c = E/((1.0+nu)*(1.0-2.0*nu))
    [ c*(1.-nu)      c*nu        c*nu             0.0             0.0             0.0
          c*nu   c*(1.-nu)       c*nu             0.0             0.0             0.0
          c*nu       c*nu    c*(1.-nu)            0.0             0.0             0.0
           0.0        0.0         0.0   c*(1.0-2.0*nu)            0.0             0.0
           0.0        0.0         0.0             0.0   c*(1.0-2.0*nu)            0.0
           0.0        0.0         0.0             0.0             0.0   c*(1.0-2.0*nu) ]
end

function mount_D(mat::ElasticSolid, ipd::ElasticSolidIpData)
    return mount_De(mat.E, mat.nu)
end

function getvals(ipd::ElasticSolidIpData)
    σ  = ipd.σ
    ε  = ipd.ε
    ndim = ipd.ndim
    sr2  = 2.0^0.5

    if ndim==2;
        return [
          :sxx => σ[1],
          :syy => σ[2],
          :szz => σ[3],
          :sxy => σ[4]/sr2,
          :exx => ε[1],
          :eyy => ε[2],
          :ezz => ε[3],
          :exy => ε[4]/sr2,
          :s_m => sum(σ[1:3])/3.0 ]
      else
        return [
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
          :s_m => sum(σ[1:3])/3.0 ]
      end
end

