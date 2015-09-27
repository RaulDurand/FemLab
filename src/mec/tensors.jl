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

