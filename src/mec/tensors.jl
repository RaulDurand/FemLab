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

# Tensor definitions using Mandel notation

export Tensor2, Tensor4

typealias Tensor2 Array{Float64,1}
typealias Tensor4 Array{Float64,2}

tensor2() = zeros(6)
tensor4() = zeros(6,6)

const tI  = [1., 1., 1., 0., 0., 0.]

import Base.trace
trace(T::Tensor2) = sum(T[1:3])
J1(T::Tensor2) = sum(T[1:3])
J2(T::Tensor2) = 0.5*dot(T,T)


const Psd = [  
    2/3. -1/3. -1/3. 0. 0. 0.
   -1/3.  2/3. -1/3. 0. 0. 0.
   -1/3. -1/3.  2/3. 0. 0. 0.
      0.    0.    0. 1. 0. 0.
      0.    0.    0. 0. 1. 0.
      0.    0.    0. 0. 0. 1. ]

const Isym = eye(6)

function dev(T::Tensor2) # deviatoric tensor
    return Psd*T
end

function J2D(T::Tensor2)
    return J2(Psd*T)
end


function J3D(T::Tensor2)
    return J3(Psd*T)
end

function principal(T::Tensor2)
    sr2 = √2.
    # full notation
    F = [ T[1]      T[4]/sr2  T[6]/sr2 ;
          T[4]/sr2  T[2]      T[5]/sr2 ;
          T[6]/sr2  T[5]/sr2  T[3]     ]
    L, V = eig(F, permute=false, scale=false)

    # force a clockwise system
    if norm( cross(V[:,2], V[:,3]) - V[:,1] ) > 1e-5
        L[2], L[3] = L[3], L[2]
        V[:,2], V[:,3] = V[:,3], V[:,2]
    end

    # find max value
    val, idx = findmax(L)
    shift = 1 - idx
    L = circshift(L, shift)
    V = circshift(V, (0,shift))
    return L, V
end


const V2M = [ 1., 1., 1., √2., √2., √2. ]
const M2V = [ 1., 1., 1., √.5, √.5, √.5 ]

function tfull(T::Tensor2)
    t1, t2, t3, t4, t5, t6 = T.*M2V
    return [ 
        t1 t4 t6
        t4 t2 t5
        t6 t5 t3 ]
end

function dyad(T1::Tensor2, T2::Tensor2)
    return T1 * T2'
end

⊗ = dyad

function inner(T1::Tensor4, T2::Tensor4)
    return sum(T1 .* T2)
end

function inner(T1::Tensor4, T2::Tensor2)
    return T1 * T2
end

function inner(T1::Tensor2, T2::Tensor4)
    return T2*T1
end

function inner(T1::Tensor2, T2::Tensor4, T3::Tensor2)
    return sum(T2*T1 .* T3)
end

∷ = inner

function rotation4(V::Array{Float64,2}, T::Tensor4)
    l1, m1, n1 = V[:,1]
    l2, m2, n2 = V[:,2]
    l3, m3, n3 = V[:,3]

    sq2 = √2.0

    T[1,1] =     l1*l1;  T[1,2] =     m1*m1;  T[1,3] =     n1*n1;  T[1,4] =   sq2*l1*m1;  T[1,5] =   sq2*m1*n1;  T[1,6] =   sq2*n1*l1     
    T[2,1] =     l2*l2;  T[2,2] =     m2*m2;  T[2,3] =     n2*n2;  T[2,4] =   sq2*l2*m2;  T[2,5] =   sq2*m2*n2;  T[2,6] =   sq2*n2*l2     
    T[3,1] =     l3*l3;  T[3,2] =     m3*m3;  T[3,3] =     n3*n3;  T[3,4] =   sq2*l3*m3;  T[3,5] =   sq2*m3*n3;  T[3,6] =   sq2*n3*l3     
    T[4,1] = sq2*l1*l2;  T[4,2] = sq2*m1*m2;  T[4,3] = sq2*n1*n2;  T[4,4] = l1*m2+l2*m1;  T[4,5] = m1*n2+m2*n1;  T[4,6] = l1*n2+l2*n1     
    T[5,1] = sq2*l2*l3;  T[5,2] = sq2*m2*m3;  T[5,3] = sq2*n2*n3;  T[5,4] = l2*m3+l3*m2;  T[5,5] = m2*n3+m3*n2;  T[5,6] = l2*n3+l3*n2     
    T[6,1] = sq2*l3*l1;  T[6,2] = sq2*m3*m1;  T[6,3] = sq2*n3*n1;  T[6,4] = l3*m1+l1*m3;  T[6,5] = m3*n1+m1*n3;  T[6,6] = l3*n1+l1*n3 
end
