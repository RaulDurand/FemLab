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

mutable struct KotsovosIpState<:IpState
    ndim::Int64
    σ::Tensor2
    ε::Tensor2
    ncracks::Int64
    D::Array{Float64,2}
    V1::Array{Float64,1} # first  crack direction
    V2::Array{Float64,1} # second crack direction
    V3::Array{Float64,1} # third  crack direction
    function KotsovosIpState(ndim=3) 
        this = new(ndim)
        this.σ   = zeros(6)
        this.ε   = zeros(6)
        this.ncracks = 0
        this.D = zeros(6,6)
        this.V1 = zeros(3)
        this.V2 = zeros(3)
        this.V3 = zeros(3)
        return this
    end
end

mutable struct Kotsovos<:AbsSolid
    E::Float64
    nu::Float64
    β::Float64
    fc::Float64
    ft::Float64
    De::Tensor4

    function Kotsovos(prms::Dict{Symbol,Float64})
        return Kotsovos(;prms...)
    end

    function Kotsovos(;E=NaN, nu=0.0, beta=0.5, fc=0.0, ft=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert 0.0<beta<=1.0
        @assert fc>=0.0  # if fc==0 => no cracks in compression
        if ft == 0.0
            ft = 0.1*fc
        end
        @assert ft>0.0

        this     = new(E, nu, beta, fc, ft)
        this.De  = calcDe(E, nu)
        return this 
    end
end

# Create a new instance of Ip data
new_ip_state(mat::Kotsovos, ndim::Int) = KotsovosIpState(ndim)

function set_state(ipd::KotsovosIpState; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig.*V2M
    else
        if length(sig)!=0; error("Kotsovos: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*V2M
    else
        if length(eps)!=0; error("Kotsovos: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::Kotsovos, ipd::KotsovosIpState)
    if ipd.ncracks == 0
        ipd.D = mat.De # defines D as elastic for this ip
    end
    return ipd.D
end

function get_basicD(mat::Kotsovos, ipd::KotsovosIpState)
    # TODO: check elastic part results with the elastic model
    ncracks = ipd.ncracks
    nu = mat.nu
    β  = mat.β
    G = 0.5*mat.E/(1+nu)
    λ = mat.E*nu/(1+nu)/(1-2*nu)
    D = zeros(6,6)

    if ncracks==0
        D = copy(mat.De)
    elseif ncracks==1
        # aligned with the x' direction
        D[2,2] = 2*G+λ
        D[3,3] = 2*G+λ
        D[2,3] = λ
        D[3,2] = λ
        D[4,4] = β*G
        D[5,5] = G
        D[6,6] = β*G
    elseif ncracks==2
        # aligned with the z'' direction
        D[3,3] = 2*G+λ
        D[4,4] = β*G
        D[5,5] = β*G
        D[6,6] = β*G
    elseif ncracks==3
        D[4,4] = β*G
        D[5,5] = β*G
        D[6,6] = β*G
    end

    # fix for Mandel notation
    D[4,4] *= 2.0
    D[5,5] *= 2.0
    D[6,6] *= 2.0

    return D
end

function stress_update(mat::Kotsovos, ipd::KotsovosIpState, Δε::Array{Float64,1})
    # initial stress
    σini = ipd.σ
    #initial number of cracks
    nini = ipd.ncracks

    # rotation matrix
    T = zeros(6,6)

    # trial stress
    ipd.D = calcD(mat, ipd)
    σtr   = ipd.σ + inner(ipd.D, Δε)

    # principal stresses
    Sig, EV = principal_dir(σtr)

    # check for cracking in traction
    #@show Sig
    if Sig[1]>mat.ft
        ipd.ncracks += 1
    end
    if Sig[2]>mat.ft
        ipd.ncracks += 1
    end
    if Sig[3]>mat.ft
        ipd.ncracks += 1
    end

    # check for cracking in compression
    if mat.fc != 0.0 && minimum(Sig) < -mat.fc
        ipd.ncracks = 3
    end

    # check if there are new cracks
    ipd.ncracks = min(ipd.ncracks, 3)
    newcracks   = nini != ipd.ncracks

    # update D matrix
    if newcracks
        if ipd.ncracks == 1
            ipd.V1 = EV[:,1]
            bD = get_basicD(mat, ipd)
            rotation4(EV, T)
            ipd.D = T'*bD*T
        end

        if ipd.ncracks == 2

            if nini == 1
                bD = get_basicD(mat, ipd)
                # align uncracked direction with the z'' direction
                ipd.V2 = EV[:,1]
                Vz = cross(ipd.V1, ipd.V2)
                Vx = cross(ipd.V2, Vz)
                Vy = cross(Vz, Vx)
                EV = hcat(Vx, Vy, Vz)

                rotation4(EV, T)
                ipd.D = T'*bD*T
            end

            if nini == 0
                bD = get_basicD(mat, ipd)

                # uncracked direction is already aligned with the z'' direction
                rotation4(EV, T)
                ipd.D = T'*bD*T
            end
        end

        if ipd.ncracks == 3
            ipd.D = get_basicD(mat, ipd)
        end
    end

    #if newcracks
        ipd.ε += Δε
        ipd.σ = inner(ipd.D, ipd.ε) # important!
        Δσ    = ipd.σ - σini
    #else
        #ipd.ε += Δε
        #ipd.σ = σtr
        #Δσ    = σtr - σini
    #end

    return Δσ 
end

function elem_vals(mat::Kotsovos, elem::Element)
    ncracks = maximum( [ ip.data.ncracks for ip in elem.ips] )

    return Dict(
     :ncracks => ncracks
    )
end

function ip_state_vals(mat::Kotsovos, ipd::KotsovosIpState)
    σ  = ipd.σ
    ε  = ipd.ε
    ndim = ipd.ndim
    sr2  = √2.

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
          :p   => sum(σ[1:3])/3.0 )
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
          :ev  => trace(ε),
          :p   => trace(σ)/3.0
          )
      end
end
