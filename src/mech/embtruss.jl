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


mutable struct EmbTruss<:AbsEmbTruss
    trussmat::Truss  # uses the material for elastic truss

    function EmbTruss(prms::Dict{Symbol,Float64})
        return EmbTruss(;prms...)
    end

    function EmbTruss(;E::Number=1.0, A::Number=1.0)
        this = new()
        this.trussmat= Truss(E=E, A=A)
        return this
    end
end

# Create a new instance of Ip data
new_ip_state(mat::EmbTruss, ndim::Int) = TrussIpState(ndim)

# Note. Not need to define function set_state because it uses Truss model with TrussIpState

function stress_update(mat::EmbTruss, ipd::TrussIpState, Δε::Float64)
    return stress_update(mat.trussmat, ipd, Δε)
end

function ip_state_vals(mat::EmbTruss, ipd::TrussIpState)
    return ip_state_vals(mat.trussmat, ipd)
end
