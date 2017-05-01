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

export Mazars
export set_state

type MazarsIpData<:IpData
    ndim::Int64
    σ::Tensor2
    ε::Tensor2
    φt::Float64
    φc::Float64
    φ::Float64  # damage
    ε̅max::Float64
    D::Array{Float64,2}
    function MazarsIpData(ndim=3) 
        this = new(ndim)
        this.σ = zeros(6)
        this.ε = zeros(6)
        this.φ = 0.0
        this.ε̅max = 0.0
        this.D = zeros(6,6)
        return this
    end
end

type Mazars<:AbsSolid
    E ::Float64
    nu::Float64
    ε̅0::Float64
    At::Float64
    Bt::Float64
    Ac::Float64
    Bc::Float64
    De::Tensor4
    invDe::Tensor4

    function Mazars(prms::Dict{Symbol,Float64})
        return Mazars(;prms...)
    end

    function Mazars(;E=NaN, nu=0.0, eps0=NaN, At=NaN, Bt=NaN, Ac=NaN, Bc=NaN)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert At>0.0
        @assert Ac>0.0
        @assert Bt>0.0
        @assert Bc>0.0

        this     = new(E, nu, eps0, At, Bt, Ac, Bc)
        this.De  = calcDe(E, nu)
        this.invDe  = inv(this.De)
        return this 
    end
end

# Create a new instance of Ip data
new_ipdata(mat::Mazars, ndim::Int) = MazarsIpData(ndim)

function set_state(ipd::MazarsIpData; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig.*V2M
    else
        if length(sig)!=0; error("Mazars: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*V2M
    else
        if length(eps)!=0; error("Mazars: Wrong size for strain array: $eps") end
    end
end


# special functions
pos(x) = (abs(x)+x)/2.0
neg(x) = (-abs(x)+x)/2.0


function calcD2(mat::Mazars, ipd::MazarsIpData)
    # There is something wrong with the derivatives here

    if ipd.φ <= 0.0
        # Elastic constitutive matrix
        return mat.De
    else
        # Principal stresses and principal directions
        σp, V = principal(ipd.σ)
        σp = [ σp; zeros(3) ]
        # Eigen vectors
        p1 = V[:, 1]
        p2 = V[:, 2]
        p3 = V[:, 3]
        P  = [ matrix2Mandel(p1*p1'), matrix2Mandel(p2*p2'), matrix2Mandel(p3*p3') ]

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)
        Dei= inv(mat.De)
        εt = Dei*σt
        εc = Dei*σc
        εv = sum(pos.(εt)) + sum(pos.(εc))

        # Equivalent strain scalar
        εp, V = principal(ipd.ε)
        ε̅ = √sum( pos(εp[i])^2 for i=1:3 )

        # Tensile and compression damage weights
        αt =  sum(pos.(εt))/εv
        αc =  sum(pos.(εc))/εv

        # Constitutive matrix calculation
        dφtdε̅ = (1.0-mat.At)*mat.ε̅0*ε̅^-2 + mat.At*mat.Bt*exp( -mat.Bt*(ε̅-mat.ε̅0) )
        dφcdε̅ = (1.0-mat.Ac)*mat.ε̅0*ε̅^-2 + mat.Ac*mat.Bc*exp( -mat.Bc*(ε̅-mat.ε̅0) )
        dε̅dε  = sum( pos(εp[i])*(1+sign(εp[i]))*P[i] for i=1:3 ) / (2*ε̅)
        #dε̅dε  = [ dε̅dε; zeros(3) ]
        dφdε  = (αt*dφtdε̅ + αc*dφcdε̅) * dε̅dε

        #@show dφdε'*mat.De
        D     = (1.0 - ipd.φ)*mat.De - (dφdε'*mat.De)'*ipd.ε'

        return D
    end
end


function calcD(mat::Mazars, ipd::MazarsIpData)
    if ipd.φ <= 0.0
        return mat.De
    else
        D = (1.0 - ipd.φ)*mat.De 
        return D
    end
end


function calc_principal(T::Array{Float64,1})
    sr2 = √2.
    # full notation
    F = [ T[1]      T[4]/sr2  T[6]/sr2 ;
          T[4]/sr2  T[2]      T[5]/sr2 ;
          T[6]/sr2  T[5]/sr2  T[3]     ]
    L, V = eig(F, permute=false, scale=false)

    return L
end


function stress_update(mat::Mazars, ipd::MazarsIpData, Δε::Array{Float64,1})
    σini  = ipd.σ
    ipd.ε = ipd.ε + Δε

    # Principal stresses tensor
    εp, V = principal(ipd.ε)

    # Equivalent strain scalar
    ε̅ = √sum( pos(εp[i])^2 for i=1:3 )
    ipd.ε̅max = max(ipd.ε̅max, mat.ε̅0)

    if ε̅ < ipd.ε̅max  # linear-elastic increment
        ipd.σ = (1.0 - ipd.φ)*mat.De*ipd.ε
    else # increment with damage
        ipd.ε̅max = ε̅

        # Principal stresses and principal directions
        #norm(ipd.σ)==0.0 && error("Mazars: At least one increment inside the elastic region is required")
        σp, V = principal(ipd.σ)
        σp = [ σp; zeros(3) ]

        # Damage calculation
        #@show mat.At
        #@show mat.Bt
        φt = 1.0 - (1-mat.At)*mat.ε̅0/ε̅ - mat.At/exp(mat.Bt*(ε̅-mat.ε̅0))
        φc = 1.0 - (1-mat.Ac)*mat.ε̅0/ε̅ - mat.Ac/exp(mat.Bc*(ε̅-mat.ε̅0))

        ipd.φt = max(ipd.φt, φt)
        ipd.φc = max(ipd.φc, φc)

        # Tensile and compression tensors
        σt = pos.(σp)
        σc = neg.(σp)

        εt = mat.invDe*σt
        εc = mat.invDe*σc
        εv = sum(pos(εt[i]) + pos(εc[i]) for i=1:3)

        # Tensile and compression damage weights
        αt =  sum(pos(εt[i]) for i=1:3)/εv
        αc =  sum(pos(εc[i]) for i=1:3)/εv
        if εv==0.0
            αt = αc = 0.5
        end

        #@show αt
        #@show αc
        #@show φt
        #@show φc

        # Damage variable
        #ipd.φ = αt*φt + αc*φc
        #φ = αt*φt + αc*φc
        φ = αt*ipd.φt + αc*ipd.φc
        ipd.φ = max(φ, ipd.φ)

        # Total stress and stress increment
        ipd.σ = (1.0 - ipd.φ)*mat.De*ipd.ε
    end

    Δσ    = ipd.σ - σini
    return Δσ 
end

function elem_vals(mat::Mazars, elem::Element)
    return Dict()
end

function getvals(mat::Mazars, ipd::MazarsIpData)
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
          :p   => trace(σ)/3.0,
          :dam => ipd.φ,
          :damt => ipd.φt,
          :damc => ipd.φc,
          :eq => ipd.ε̅max
          )
      end
end
