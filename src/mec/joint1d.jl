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


export Joint1D
abstract AbsJoint1D<:Mechanical

type Joint1DIpData<:IpData
    ndim::Int
    sig ::Array{Float64,1}
    eps ::Array{Float64,1}
    function Joint1DIpData(ndim=3)
        this = new(ndim)
        this.sig = zeros(3)
        this.eps = zeros(3)
        return this
    end
end

type Joint1D<:AbsJoint1D
    ks::Float64
    kn::Float64
    h ::Float64    # section perimeter
    new_ipdata::DataType

    function Joint1D(prms::Dict{Symbol,Float64})
        return  Joint1D(;prms...)
    end

    function Joint1D(;ks=NaN, kn=NaN, h=NaN, A=NaN, dm=NaN)
        # A : section area
        # dm: section diameter
        # h : section perimeter
        @assert ks>=0
        @assert kn>=0
        @assert (h>0 || A>0 || dm>0)

        if isnan(h) 
            if A>0
                h = 2.0*(A*pi)^0.5
            else
                h = pi*dm
            end
        end
        @assert h>0

        this = new(ks, kn, h)
        this.new_ipdata = Joint1DIpData
        return this
    end
end

function set_state(ipd::Joint1DIpData, sig=zeros(0), eps=zeros(0))
    if length(sig)==3
        ipd.sig[:] = sig
    else
        if length(sig)!=0; error("MecElasticSolid: Wrong size for stress array: $sig") end
    end
    if length(eps)==3
        ipd.eps[:] = eps
    else
        if length(eps)!=0; error("MecElasticSolid: Wrong size for strain array: $eps") end
    end
end

function calcD(mat::Joint1D, ipd::Joint1DIpData)
    ks = mat.ks
    kn = mat.kn
    if ipd.ndim==2
        return [  ks  0.0 
                 0.0   kn ]
    else
        return  [  ks  0.0  0.0
                  0.0   kn  0.0
                  0.0  0.0   kn ]
    end
end


function stress_update(mat::Joint1D, ipd::Joint1DIpData, deps)
    D = calcD(mat, ipd)
    dsig     = D*deps

    ipd.eps[1:ipd.ndim] += deps
    ipd.sig[1:ipd.ndim] += dsig
    return dsig
end


function mountB(mat::AbsJoint1D, elem::Element, R, Ch, Ct, B)
    # Calculates the matrix that relates nodal displacements with relative displacements
    # ==================================================================================

    # B = T* [NN*MM  -NN]      ndim x ndim*(m+n)

    # where
    # T is a direction cosines matrix
    # NN is a matrix containing truss element shape functions
    # evaluated at the point of interest R.
    # MM is a matrix containing tresspased element shape functions
    # evaluated at the n truss nodal points.

    #          [ M_11*I M_21*I ... M_m1*I]
    #     MM = [ M_12*I M_22*I ... M_m2*I]
    #          [ M_13*I M_23*I ... M_mn*I] ndim*n x ndim*m

    # where
    # M_12 is the first shape function from tresspased element
    # evaluated at the second node of truss element.
    # I is a ndim x ndim identity matrix


    ndim = elem.ndim
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    nnodes  = length(elem.nodes)
    ntnodes = length(bar.nodes)
    nbnodes = length(hook.nodes)
    D = deriv_func(bar.shape, R)
    J = D*Ct
    T = mount_T(J)

    # Mount NN matrix
    N = shape_func(bar.shape, R)
    NN = hcat([ Ni*eye(ndim) for Ni in N  ]...)

    # Mount MM matrix
    stack = Array{Float64,2}[]
    for i=1:ntnodes
        Xj = bar.nodes[i].X
        R  = inverse_map(hook.shape, Ch, Xj)
        M  = shape_func(hook.shape, R)
        for Mi in M
            push!(stack, Mi*eye(ndim))
        end
    end
    MM = hvcat(nbnodes, stack...)

    B[:] = -T*[ NN*MM  -NN ]
    detJ = norm(J)
    return detJ
end

function elem_jacobian(mat::AbsJoint1D, elem::Element)
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    Ch = getcoords(hook)
    Ct = getcoords(bar)

    K  = zeros(nnodes*ndim, nnodes*ndim)
    B  = zeros(ndim, nnodes*ndim)
    DB = zeros(ndim, nnodes*ndim)

    for ip in elem.ips
        detJ = mountB(elem.mat, elem, ip.R, Ch, Ct, B) #TODO: include B matrix construction here
        D    = calcD(elem.mat, ip.data)
        D[1,1]*=mat.h
        coef = detJ*ip.w
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end
    return K
end

function update!(mat::AbsJoint1D, elem::Element, DU::Array{Float64,1}, DF::Array{Float64,1})
    ndim   = elem.ndim
    nnodes = length(elem.nodes)
    mat    = elem.mat
    map    = get_map(elem)

    dF = zeros(nnodes*ndim)
    dU = DU[map]
    B  = zeros(ndim, nnodes*ndim)

    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]
    Ct   = getcoords(bar)
    Ch   = getcoords(hook)
    for ip in elem.ips
        detJ = mountB(elem.mat, elem, ip.R, Ch, Ct, B)
        D    = calcD(mat, ip.data)
        deps = B*dU
        dsig = stress_update(mat, ip.data, deps)
        coef = detJ*ip.w
        dsig[1]  *= mat.h
        @gemv dF += coef*B'*dsig
    end

    # Update global vector
    DF[map] += dF
end

#function getvals(ipd::Joint1DIpData)
    #return Dict(
      #:ur   => ipd.eps[1] ,
      #:tau  => ipd.sig[1] )
#end

function getvals(mat::Joint1D, ipd::Joint1DIpData)
    return Dict(
      :ur   => ipd.eps[1] ,
      :tau  => ipd.sig[1] )
end

function node_and_elem_vals(mat::AbsJoint1D, elem::Element)
    ndim = elem.ndim
    node_vals = Dict{Symbol, Array{Float64,1}}()
    elem_vals = Dict{Symbol, Float64}()

    for key in (:ux, :uy, :uz)[1:ndim]
        node_vals[key] = [node.dofdict[key].U for node in elem.nodes]
    end

    # Elem vals
    all_ip_vals = [ getvals(mat, ip.data) for ip in elem.ips ]
    labels      = keys(all_ip_vals[1])
    nips        = length(elem.ips) 

    # matrix with all ip values (nip x nvals)
    IP = vcat([ [values(all_ip_vals[i])...]' for i=1:nips]...)

    hook = elem.linked_elems[1]
    bar  = elem.linked_elems[2]

    E = extrapolator(bar.shape, nips)
    N = E*IP # (nnodes x nvals)

    nhnodes = length(hook.nodes)
    N = vcat( zeros(nhnodes, length(labels)), N )

    # Filling nodal and elem vals
    for (i,key) in enumerate(labels)
        node_vals[key] = N[:,i]
    end

    return node_vals, elem_vals

end

