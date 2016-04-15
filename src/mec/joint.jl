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


export Joint
abstract AbsJoint<:Mechanical

type JointIpData<:IpData
    ndim::Int
    sig ::Array{Float64,1}
    eps ::Array{Float64,1}
    function JointIpData(ndim=3)
        this = new(ndim)
        this.sig = zeros(3)
        this.eps = zeros(3)
        return this
    end
end

type Joint<:AbsJoint
    ks::Float64
    kn::Float64
    new_ipdata::DataType

    function Joint(prms::Dict{Symbol,Float64})
        return  Joint(;prms...)
    end

    function Joint(;ks=NaN, kn=NaN)
        @assert ks>=0
        @assert kn>=0

        this = new(ks, kn)
        this.new_ipdata = JointIpData
        return this
    end
end

function set_state(ipd::JointIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.sig[:] = sig
    else
        if length(sig)!=0; error("Joint: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.eps[:] = eps
    else
        if length(eps)!=0; error("Joint: Wrong size for strain array: $eps") end
    end
end

function mountD(mat::Joint, ipd::JointIpData)
    ks = mat.ks
    kn = mat.kn
    if ipd.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   ks  0.0
                  0.0  0.0   kn ]
    end
end

function stress_update(mat::Joint, ipd::JointIpData, Δu)
    D  = mountD(mat, ipd)
    Δσ = D*Δu

    ipd.eps[1:ipd.ndim] += Δu
    ipd.sig[1:ipd.ndim] += Δσ
    return Δσ
end

function matrixT(J::Matrix{Float64})
    if size(J,1)==2
        L1 = vec(J[1,:])
        L2 = vec(J[2,:])
        L3 = cross(L1, L2)
        L1 /= norm(L1)
        L2 /= norm(L2)
        L3 /= norm(L3)
        return [L1 L2 L3]'
    else
        L1 /= vec(J)/norm(J)
        L2  = [ -L1[2],  L1[1] ]
        return [L1 L2]'
    end
end

function elem_jacobian(mat::AbsJoint, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    bshape = get_basic_shape(elem.shape)

    C = getcoords(elem)
    B = zeros(ndim, nnodes*ndim)
    K = zeros(nnodes*ndim, nnodes*ndim)

    DB = zeros(ndim, nnodes*ndim)
    J  = zeros(ndim-1, ndim)
    NN = zeros(ndim, nnodes*ndim)

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = deriv_func(bshape, ip.R)
        N    = shape_func(bshape, ip.R)
        @gemm J = dNdR*C
        detJ = norm2(J)

        # compute B matrix
        T    = matrixT(J)
        NN[:,:] = 0.0  # NN = [ -N[]  N[] ]
        for i=1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end
        @gemm B = T*NN

        # compute K
        coef = detJ*ip.w
        D    = mountD(elem.mat, ip.data)
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end

    return K
end

function update!(mat::AbsJoint, elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    bshape = get_basic_shape(elem.shape)
    mat    = elem.mat
    map    = get_map(elem)

    dF = zeros(nnodes*ndim)
    dU = DU[map]
    C = getcoords(elem)
    B = zeros(ndim, nnodes*ndim)

    DB = zeros(ndim, nnodes*ndim)
    J  = zeros(ndim-1, ndim)
    NN = zeros(ndim, nnodes*ndim)
    Δu = zeros(3)

    C    = getcoords(elem)

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = deriv_func(bshape, ip.R)
        N    = shape_func(bshape, ip.R)
        @gemm J = dNdR*C
        detJ = norm2(J)

        # compute B matrix
        T    = matrixT(J)
        NN[:,:] = 0.0  # NN = [ -N[]  N[] ]
        for i=1:hnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end
        @gemm B = T*NN

        # internal force
        @gemv Δu = B*dU
        Δσ   = stress_update(mat, ip.data, Δu)
        coef = detJ*ip.w
        @gemv dF += coef*B'*Δσ
    end

    # Update global vector
    DF[map] += dF
end

function getvals(mat::Joint, ipd::JointIpData)
    if ipd.ndim == 2
        return Dict(
          :s11  => ipd.sig[1] ,
          :s12  => ipd.sig[2] )
    else
        return Dict(
          :s11  => ipd.sig[1] ,
          :s12  => ipd.sig[2] ,
          :s13  => ipd.sig[3] )
    end
end

function node_and_elem_vals(mat::AbsJoint, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    # Elem vals
    #all_ip_vals = [ getvals(mat, ip.data) for ip in elem.ips ]
    #labels      = keys(all_ip_vals[1])
    #nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    #IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    #E = extrapolator(get_basic_shape(shape), nips)
    #N = E*IP # (nnodes x nvals)

    # Filling nodal and elem vals
    #for (i,key) in enumerate(labels)
        #node_vals[key] = N[:,i]
    #end

    return node_vals, elem_vals

end

