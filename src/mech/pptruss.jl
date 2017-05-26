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

type PPTrussIpData<:IpData
    ndim::Int
    σ::Float64
    ε::Float64
    εpa::Float64
    Δγ ::Float64

    function PPTrussIpData(ndim=3)
        this = new(ndim)
        this.σ = 0.0
        this.ε = 0.0
        this.εpa = 0.0
        this.Δγ  = 0.0
        return this
    end
end

type PPTruss<:AbsTruss
    E::Float64
    A::Float64
    σy0::Float64
    H::Float64

    function PPTruss(prms::Dict{Symbol,Float64})
        return  PPTruss(;prms...)
    end

    function PPTruss(;E=NaN, A=NaN, sig_y=NaN, H=0.0)
        @assert E>0
        @assert A>0
        @assert sig_y>0
        this = new(E, A, sig_y, H)
        this
    end
end

# Create a new instance of Ip data
new_ipdata(mat::PPTruss, ndim::Int) = PPTrussIpData(ndim)

function set_state(ipd::PPTrussIpData, σ=NaN, ε=NaN)
    if !isnan(σ); ipd.σ = σ end
    if !isnan(ε); ipd.ε = ε end
end

function yield_func(mat::PPTruss, ipd::PPTrussIpData, σ::Float64)
    σya = mat.σy0 + mat.H*ipd.εpa
    return abs(σ) - σya
end

function calcD(mat::PPTruss, ipd::PPTrussIpData)
    if ipd.Δγ == 0.0
        return mat.E
    else
        E, H = mat.E, mat.H
        return E*H/(E+H)
    end
end

function stress_update(mat::PPTruss, ipd::PPTrussIpData, Δε::Float64)
    E, H = mat.E, mat.H
    σini = ipd.σ
    σtr    = σini + E*Δε
    ftr    = yield_func(mat, ipd, σtr)
    ipd.Δγ = ftr>0.0? ftr/(E+H) : 0.0
    Δεp    = ipd.Δγ*sign(σtr)
    ipd.εpa += ipd.Δγ
    ipd.σ  = σtr - E*Δεp
    Δσ     = ipd.σ - σini
    ipd.ε += Δε
    return Δσ
end

#=
function getvals(ipd::PPTrussIpData)
    return Dict(
      :sa => ipd.σ,
      :ea => ipd.ε,
      :epa_ppt => ipd.εpa )
      #:Fa => ipd.σ*mat.A,
      #:A  => mat.A ]
end
=#

function getvals(mat::PPTruss, ipd::PPTrussIpData)
    return Dict(
      :sa => ipd.σ,
      :ea => ipd.ε,
      :epa_ppt => ipd.εpa,
      :Fa => ipd.σ*mat.A,
      :A  => mat.A )
end

