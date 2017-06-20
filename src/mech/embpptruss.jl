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


mutable struct EmbPPTruss<:AbsEmbTruss
    trussmat::PPTruss  # uses the material for elastic truss

    function EmbPPTruss(prms::Dict{Symbol,Float64})
        return EmbPPTruss(;prms...)
    end

    function EmbPPTruss(;E=NaN, A=NaN, sig_y=NaN, H=0.0)
        this = new()
        this.trussmat = PPTruss(E=E, A=A, sig_y=sig_y, H=H)
        return this
    end
end

# Create a new instance of Ip data
new_ip_state(mat::EmbPPTruss, ndim::Int) = PPTrussIpState(ndim)

# Obs. Not need to define function set_state because it uses Truss model with PPTrussIpState

function stress_update(mat::EmbPPTruss, ipd::PPTrussIpState, Δε::Float64)
    return stress_update(mat.trussmat, ipd, Δε)
end

function ip_state_vals(mat::EmbPPTruss, ipd::PPTrussIpState)
    return ip_state_vals(mat.trussmat, ipd)
end

function calcD(mat::EmbPPTruss, ips::PPTrussIpState)
    return mat.trussmat.E
end


