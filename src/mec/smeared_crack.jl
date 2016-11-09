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

export SmearedCrack
export set_state

type SmearedCrackIpData<:IpData
    ndim::Int64
    σ::Tensor2
    ε::Tensor2
    εcr::Array{Float64,1}
    ncracks::Int64
    V::Array{Float64,2} # crack direction
    function SmearedCrackIpData(ndim=3) 
        this = new(ndim)
        this.σ   = zeros(6)
        this.ε   = zeros(6)
        this.εcr = zeros(3)
        this.ncracks = 0
        return this
    end
end

type SmearedCrack<:Mechanical
    E ::Float64  # concrete Young's modulus
    nu::Float64  # Poisson ratio
    fc::Float64  # compressive strenght of concrete
    ft::Float64  # tensile strenght
    Gf::Float64  # Fracture energy
    ξ1::Float64  # horizontal normalized coordinate for a point in σn vs εn graph
    ξ2::Float64
    α1::Float64  # vertical normalized coordinate for a point in σn vs εn graph
    α2::Float64
    p1::Float64  # exponent in β formula
    h ::Float64  # band-width of a crack
    De::Tensor4  # Elastic tensor
    new_ipdata::DataType #TODO: use an instance of SmearedCrackIpData instead, later use copies of it

    function SmearedCrack(prms::Dict{Symbol,Float64})
        return SmearedCrack(;prms...)
    end

    function SmearedCrack(;E=NaN, nu=0.0, fc=0.0, ft=0.0, Gf=0.0, xi1=0.0, xi2=0.0, al1=0.0, al2=0.0, p1=0.0, h=0.0)
        @assert E>0.0
        @assert 0.0<=nu<0.5
        @assert fc>=0.0  # if fc==0 => no cracks in compression
        if ft == 0.0
            ft = 0.1*fc
        end
        @assert ft>0.0
        @assert Gf>0.0

        #@assert h <= Gf*E/(b*ft^2)

        this     = new(E, nu, fc, ft, Gf, xi1, xi2, al1, al2, p1, h)
        this.new_ipdata = SmearedCrackIpData
        this.De  = zeros(6,6)
        setDe(E, nu, this.De) # elastic tensor

        α1 = al1
        α2 = al2
        ξ1 = xi1
        ξ2 = xi2

        k1 = (1-α1)*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*ξ1)
        k2 = (α1-α2)*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*(ξ2-ξ1))
        k3 = α2*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*(1-ξ2))
        k4 = 2./(ξ1 + α1*ξ2 - α2*ξ1 + α2)
        εnul = k4*Gf/(ft*h)

        b = max(k1,k2,k3,k4)
        @assert h<=Gf*E/(b*ft^2)

        return this 
    end
end

function set_state(ipd::SmearedCrackIpData; sig=zeros(0), eps=zeros(0))
    if length(sig)==6
        ipd.σ[:] = sig.*V2M
    else
        if length(sig)!=0; error("SmearedCrack: Wrong size for stress array: $sig") end
    end
    if length(eps)==6
        ipd.ε[:] = eps.*V2M
    else
        if length(eps)!=0; error("SmearedCrack: Wrong size for strain array: $eps") end
    end
end

function calcT(V::Array{Float64,2})
    l1, m1, n1 = V[:,1]
    l2, m2, n2 = V[:,2]
    l3, m3, n3 = V[:,3]

    sq2 = √2.0 # form Mandel notation
    T = zeros(3,6)

    T[1,1] =     l1*l1;  T[1,2] =     m1*m1;  T[1,3] =     n1*n1;  T[1,4] =   sq2*l1*m1;  T[1,5] =   sq2*m1*n1;  T[1,6] =   sq2*n1*l1     
    T[2,1] = sq2*l1*l2;  T[2,2] = sq2*m1*m2;  T[2,3] = sq2*n1*n2;  T[2,4] = l1*m2+l2*m1;  T[2,5] = m1*n2+m2*n1;  T[2,6] = l1*n2+l2*n1     
    T[3,1] = sq2*l3*l1;  T[3,2] = sq2*m3*m1;  T[3,3] = sq2*n3*n1;  T[3,4] = l3*m1+l1*m3;  T[3,5] = m3*n1+m1*n3;  T[3,6] = l3*n1+l1*n3 
    return T

end


function calcσn(mat::SmearedCrack, ipd::SmearedCrackIpData, εn)
    if εn<0
        return 0.0
    end

    T   = calcT(ipd.V)

    α1, α2, ξ1, ξ2 = mat.α1, mat.α2, mat.ξ1, mat.ξ2

    k4 = 2./(ξ1 + α1*ξ2 - α2*ξ1 + α2)
    εnul = k4*mat.Gf/(mat.ft*mat.h)

    k=0

    if εn <= ξ1*εnul # k1
        k = (1-α1)*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*ξ1)
    elseif εn <= ξ2*εnul # k2
        k = (α1-α2)*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*(ξ2-ξ1))
    elseif εn <= εnul # k2
        k = α2*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*(1-ξ2))
    elseif εn > εnul
        k = 0.0
    end
    DI = -k*mat.h*mat.ft^2/mat.Gf

    if εn <= ξ1*εnul 
        return mat.ft + DI*εn
    elseif εn <= ξ2*εnul
        return α1*mat.ft + DI*(εn-ξ1*εnul)
    elseif εn <= εnul 
        return α2*mat.ft + DI*(εn-ξ2*εnul)
    elseif εn > εnul
        return 0.0
    end

end


function calcDcr(mat::SmearedCrack, ipd::SmearedCrackIpData, εn)
    T   = calcT(ipd.V)
    #εcr = T*ipd.ε

    #εn  = ipd.εcr[1]

    #if εn<0.0
        #return zeros(3,3)
    #end

    #@show ipd.εcr

    α1, α2, ξ1, ξ2 = mat.α1, mat.α2, mat.ξ1, mat.ξ2

    k4 = 2./(ξ1 + α1*ξ2 - α2*ξ1 + α2)
    εnul = k4*mat.Gf/(mat.ft*mat.h)
    #@show εn
    #@show εnul

    k=0

    if εn <= ξ1*εnul # k1
        k = (1-α1)*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*ξ1)
    elseif εn <= ξ2*εnul # k2
        k = (α1-α2)*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*(ξ2-ξ1))
    elseif εn <= εnul # k2
        k = α2*(ξ1 + α1*ξ2 - α2*ξ1 + α2)/(2*(1-ξ2))
    elseif εn > εnul
        k = 0.0
    end
    DI = -k*mat.h*mat.ft^2/mat.Gf

    #@assert εn>0
    #if εn <= 0.0
        #@show εn
    #end

    if εn>εnul   # TODO: Check
        DII = 0.0
    else
        β   = (1-εn/εnul)^mat.p1
        if εn==0
            β=0.9999999
        end
        Gc  = mat.E/(2*(1+mat.nu))
        DII = Gc*β/(1-β)
    end

    Dcr = [ DI 0. 0.; 0. DII 0.; 0. 0. DII]
    #@show Dcr
    #exit()
    return Dcr
end

function calcD(mat::SmearedCrack, ipd::SmearedCrackIpData)
    if ipd.εcr[1] <= 0.0
        return mat.De # defines D as elastic for this ip
    else
        #@show ipd.ncracks
        T = calcT(ipd.V)
        Dcr = calcDcr(mat, ipd, ipd.εcr[1])
        Dco = mat.De
        Dcrco = Dco - Dco*T'*inv(Dcr + T*Dco*T')*T*Dco
        return Dcrco
    end
end

function stress_update(mat::SmearedCrack, ipd::SmearedCrackIpData, Δε::Array{Float64,1})
    # initial stress
    σini = ipd.σ

    #initial number of cracks
    σtr = ipd.σ + mat.De*Δε

    # principal stresses
    if ipd.ncracks == 0
        Sig, EV = principal(σtr)
        #@show Sig
        #@show EV
        #exit()
        if Sig[1] > mat.ft
            #@show σtr
            #@show Sig
            ipd.ncracks = 1
            ipd.V = EV
            #@show ipd.V
            #TT = zeros(6,6)
            #rotation4(ipd.V,TT)
            #@showm TT*σtr
            #exit()
        end
    end

    if ipd.ncracks == 0
        Δσ = σtr - σini
    else
        #@show ipd.ncracks
        T = calcT(ipd.V)
        Δεcr0 = zeros(3) 
        Δεcr0 = T*Δε*0.5

        Δεcr = zeros(3)
        Dco  = mat.De

        # Fixed point iterations to find Δεcr
        maxit = 20
        i = 0
        #println()
        εn = ipd.εcr[1]
        Dcr = calcDcr(mat, ipd, εn)
        #@showm Dcr
        Δεcr = inv(Dcr + T*Dco*T')*T*Dco*Δε

        #εn = (ipd.εcr + Δεcr)[1]
        #@show εn
        @show Δεcr
        #DII = Dcr[2,2]
        #@show DII
        if false
        for i=1:maxit
            #@show i
            εn = (ipd.εcr + Δεcr0)[1]
            #σn = calcσn(mat, ipd, εn)
            #Δσn = σn - (T*σini)[1]

            #Δτ2 = DII*Δεcr0[2]
            #Δτ3 = DII*Δεcr0[3]
            #Δσ = [Δσn, Δτ2, Δτ3]
            Δσ = Dcr*Δεcr0
            
            #@show Δσ

            #@show εn
            #@show σn
            #@show Dcr
            #if DII>0
                #s = √2.0/2
                #ipd.V = [1 0 0; 0 -s -s; 0 -s s]
                #@showm ipd.V
                #@showm Δε
                #@showm Dco*Δε
                #@showm T*Dco*Δε
                #TT = zeros(6,6)
                #rotation4(ipd.V,TT)
                #@showm TT*Dco*Δε
                #exit()
            #end
            Δεcr = inv(T*Dco*T')*( T*Dco*Δε - Δσ )
            #@show Δεcr
            #@show norm(Δεcr - Δεcr0)
            #@show norm(Δεcr - Δεcr0)
            if norm(Δεcr - Δεcr0)< 1e-6 break end
            Δεcr0 .= Δεcr  # TODO: check

            #@show Δεcr0
        end
        end

        if i==maxit
            error("SmearedCrack: Not converged")
        end

        Δσ = Dco*(Δε - T'*Δεcr)
        #Δσ = Dco*Δε - Dco*T'*Δεcr
        ipd.εcr += Δεcr
    end

    ipd.ε += Δε
    ipd.σ += Δσ

    return Δσ 
end

function elem_vals(mat::SmearedCrack, elem::Element)
    ncracks = maximum( [ ip.data.ncracks for ip in elem.ips] )

    return Dict(
     :ncracks => ncracks
    )
end

function getvals(mat::SmearedCrack, ipd::SmearedCrackIpData)
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
