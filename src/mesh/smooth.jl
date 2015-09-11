
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

# This file contains the code for smoothing meshes

export smooth!

# Returns a matrix with the cell coordinates
function cellcoords(c::Cell)
    ndim = get_ndim(c.shape)
    C = Array(Float64, length(c.points), ndim)
    for (i,point) in enumerate(c.points)
        C[i,1] = point.x
        C[i,2] = point.y
        if ndim>2
            C[i,3] = point.z
        end
    end
    return C
end

# Basic coordinates are defined considering an area/volume equal to 1.0
function basic_coords(shape::ShapeType) #check
    if shape == TRI3
        #return √2*[0. 0; 1 0 ; 0 1]
        a = 2./3.^0.25
        return [ 0. 0.; a 0.; a/2. a/2*√3 ]
    end
    if shape == TRI6
        #return √2*[0. 0; 1 0 ; 0 1; 0.5 0; 0.5 0.5; 0 0.5]
        a = 2./3.^0.25
        h = a/2*√3
        return [ 0. 0; a 0; a/2. a/2*√3; a/2. 0; 0.75*a h/2.; 0.25*a h/2. ]
    end
    if shape == QUAD4
        return [ 0. 0.; 1. 0.; 1. 1. ; 0. 1. ]
    end
    if shape == QUAD8 ; 
        return [ 0. 0; 1 0; 1 1 ; 0 1; 0.5 0; 1 0.5; 0.5 1; 0 0.5 ]
    end
    if shape == HEX8  
        return [ 0. 0. 0.; 1. 0. 0.; 1. 1. 0.; 0. 1. 0.; 0. 0. 1.; 1. 0. 1.; 1. 1. 1.; 0. 1. 1. ]
    end
    if shape == HEX20
        return [ 0.0 0.0 0.0; 
                 1.0 0.0 0.0;
                 1.0 1.0 0.0;
                 0.0 1.0 0.0;
                 0.0 0.0 1.0;
                 1.0 0.0 1.0;
                 1.0 1.0 1.0;
                 0.0 1.0 1.0;
                 0.5 0.0 0.0;
                 1.0 0.5 0.0;
                 0.5 1.0 0.0;
                 0.0 0.5 0.0;
                 0.5 0.0 1.0;
                 1.0 0.5 1.0;
                 0.5 1.0 1.0;
                 0.0 0.5 1.0;
                 0.0 0.0 0.5;
                 1.0 0.0 0.5;
                 1.0 1.0 0.5;
                 0.0 1.0 0.5 ]
    end


    error("No basic coordinates for shape $shape")
end

# Returns a rotation matrix for a cell based in their first points
function cell_orientation(cell::Cell)
    shape = cell.shape

    if shape in [ TRI3, TRI6, TRI9, TRI10, QUAD4, QUAD8, QUAD12, QUAD16 ] 
        p1 = cell.points[1]
        p2 = cell.points[2]

        l1 = p2.x - p1.x
        m1 = p2.y - p1.y
        T1 = [l1, m1]
        l1, m1 = T1/norm(T1)
        return [ l1  -m1; m1  l1 ]
    end

    if shape in [ TET4, TET10, HEX8, HEX20 ] 
        p1 = cell.points[1]
        p2 = cell.points[2]
        if shape in [ TET4, TET10 ]
            p3 = cell.points[3]
        else
            p3 = cell.points[4]
        end

        T1 = [ p2.x - p1.x, p2.y - p1.y, p2.z - p1.z ]
        T2 = [ p3.x - p1.x, p3.y - p1.y, p3.z - p1.z ]
        T3 = cross(T1, T2)
        T2 = cross(T3, T1) # redefine T2

        T1 = T1/norm(T1)
        T2 = T2/norm(T2)
        T3 = T3/norm(T3)

        return [T1 T2 T3]

    end
end


# Matrix D for the simplified FEM analysis
function matrixD(E::Float64, nu::Float64)
    c = E/((1.0+nu)*(1.0-2.0*nu))
    [ c*(1.-nu)      c*nu        c*nu             0.0             0.0             0.0
          c*nu   c*(1.-nu)       c*nu             0.0             0.0             0.0
          c*nu       c*nu    c*(1.-nu)            0.0             0.0             0.0
           0.0        0.0         0.0   c*(1.0-2.0*nu)            0.0             0.0
           0.0        0.0         0.0             0.0   c*(1.0-2.0*nu)            0.0
           0.0        0.0         0.0             0.0             0.0   c*(1.0-2.0*nu) ]
end


#include("../tools/linalg.jl")

# Matrix B for the simplified FEM analysis
function matrixB(ndim::Int, dNdX::Matx, detJ::Float64, B::Matx)
    nnodes = size(dNdX,2)
    sqr2 = √2.0
    B[:] = 0.0
    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[1,i]
            B[2,2+j*ndim] = dNdX[2,i]
            B[4,1+j*ndim] = dNdX[2,i]/sqr2; B[4,2+j*ndim] = dNdX[1,i]/sqr2
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[1,i]
            dNdy = dNdX[2,i]
            dNdz = dNdX[3,i]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,1+j*ndim] = dNdy/sqr2;   B[4,2+j*ndim] = dNdx/sqr2
            B[5,2+j*ndim] = dNdz/sqr2;   B[5,3+j*ndim] = dNdy/sqr2
            B[6,1+j*ndim] = dNdz/sqr2;   B[6,3+j*ndim] = dNdx/sqr2
        end
    end

    return detJ
end

function matrixK(cell::Cell, ndim::Int64, E::Float64, nu::Float64)
    nnodes = length(cell.points)

    C = cellcoords(cell)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)

    DB = Array(Float64, 6, nnodes*ndim)
    J = Array(Float64, ndim, ndim)
    dNdX = Array(Float64, ndim, nnodes)

    IP = get_ip_coords(cell.shape)
    D = matrixD(E, nu)

    for i=1:size(IP,1)
        R    = vec(IP[i,1:3])
        w    = IP[i,4]

        # compute B matrix
        dNdR = deriv_func(cell.shape, R)
        @gemm J = dNdR*C
        @gemm dNdX = inv(J)*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative jacobian determinant in cell $(cell.id)")
        matrixB(ndim, dNdX, detJ, B)

        # compute K
        coef = detJ*w
        @gemm DB = D*B
        @gemm K += coef*B'*DB
    end
    return K
end

function matrixK2(cell::Cell, ndim::Int64, E::Float64, nu::Float64)
    nnodes = length(cell.points)

    C = cellcoords(cell)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(6, nnodes*ndim)
    DB = Array(Float64, 6, nnodes*ndim)

    IP = get_ip_coords(cell.shape)

    D = matrixD(E, nu)
    for i=1:size(IP,1)
        R    = vec(IP[i,1:3])
        w    = IP[i,4]
        detJ = matrixB(cell, ndim, R, C, B)

        @gemm DB = D*B
        coef = detJ*w
        @gemm K += coef*B'*DB
        #K   += B'*D*B*detJ*w
    end
    return K
end

function get_map(c::Cell)
    ndim = get_ndim(c.shape)
    
    map = Int[]
    for p in c.points
        for i=1:ndim
            push!(map, (p.id-1)*ndim + i)
        end
    end

    return map
end

# Mount global stiffness matrix
function mountKg(mesh::Mesh, E::Float64, nu::Float64, A::Array{Float64,2})
    ndim = mesh.ndim
    ndof = ndim*length(mesh.points)
    R, C, V = Int64[], Int64[], Float64[]

    for c in mesh.cells
        Ke  = matrixK(c, ndim, E, nu)
        map = get_map(c)
        nr, nc = size(Ke)
        for i=1:nr
            for j=1:nc
                push!(R, map[i])
                push!(C, map[j])
                push!(V, Ke[i,j])
            end
        end
    end

    # mount A and A'
    nbc = size(A,1)
    for i = 1:nbc
        for j = 1:ndof
            val = A[i,j]
            if val==0.0
                continue
            end
            push!(R, ndof+i)
            push!(C, j)
            push!(V, val)
            push!(R, j)
            push!(C, ndof+i)
            push!(V, val)
        end
    end
    
    return sparse(R, C, V, ndof+nbc, ndof+nbc)
end

# Check if an array of faces are coplanar
# Returns the normal unitary vector of a face
# If the face is not flat returns [0,0,0]
function normal_to_faces(faces::Array{Cell, 1})
    ndim = 1 + get_ndim(faces[1].shape)
    points = Array(Point, 0)

    for f in faces
        for p in f.points
            push!(points, p)
        end
    end

    # mount coordinates matrix
    C = Array(Float64, length(points), ndim)
    for (i,p) in enumerate(points)
        C[i,1] = p.x
        C[i,2] = p.y
        if ndim>2
            C[i,3] = p.z
        end
    end

    C += 1.e-5
    I = ones(size(C,1))
    N = pinv(C)*I # best fit normal

    if norm(C*N - I)<1e-5
        return N/norm(N)
    end

    return zeros(ndim)
end

# Auxiliary structure
type sNode
    point::Point
    faces::Array{Cell}
    normal::Array{Float64}
end

# Mount global matrix A
function mountA(mesh::Mesh, fixed::Bool)
    # get border faces
    ndim = mesh.ndim
    border_fcs = get_surface(mesh.cells)

    # get all points from surface including a list of faces per point
    nodesd = Dict{Uint64, sNode}()
    for cell in border_fcs
        for point in cell.points
            hk = hash(point)
            if !haskey(nodesd, hk)
                nodesd[hk] = sNode(point, Cell[cell], Float64[])
            else
                push!(nodesd[hk].faces, cell)
            end
        end
    end
    border_nodes = [ values(nodesd)... ]

    # calc normals
    for node in border_nodes
        node.normal = normal_to_faces(node.faces)
    end

    # find the number of bcs
    n = 0  # number of bcs 
    for node in border_nodes
        if node.normal == zeros(ndim) # no coplanar faces
            n +=ndim
        else
            n +=1
        end
    end


    #  for fixed boundary
    if fixed
        n = length(border_nodes)
        A = zeros(n*ndim, length(mesh.points)*ndim)

        for i=1:n
            for j=1:ndim 
                A[ (i-1)*ndim+j, (border_nodes[i].point.id-1)*ndim+j ] = 1.0
            end
        end
    else
        # mount matrix A (Lagrange multipliers) according to bcs
        A = zeros(n, length(mesh.points)*ndim)

        baserow = 0
        for node in border_nodes
            basecol = (node.point.id-1)*ndim
            if node.normal == zeros(ndim) # no coplanar faces
                for j=1:ndim 
                    A[ baserow+j, basecol+j ] = 1.0
                end
                baserow += ndim
            else
                for j=1:ndim 
                    A[ baserow+1, basecol+j ] = node.normal[j]
                end
                baserow += 1
            end
        end
    end

    return A
end

using Base.Test

function rigid_transform(source::Array{Float64,2}, target::Array{Float64,2})
    A, B = copy(source), copy(target)
    #@test size(A) == size(B)

    n = size(A,1)
    ndim = size(A,2)

    # Centralizing both sets of points
    cA = mean(A,1)
    cB = mean(B,1)
    for i = 1:n
        for j = 1:ndim
            A[i,j] -= cA[1,j]
            B[i,j] -= cB[1,j]
        end
    end

    # Singular value decomposition
    U, S, V = svd(A'*B)

    # Rotation matrix
    R = V*U'
    #R = U*V'

    # special reflection case
    if det(R) < 0.0
       #println("Reflection detected")
       R[:2] *= -1
    end

    D = cB - cA*R'

    return R, D

end


# Mount a vector with nodal forces
function force_bc(mesh::Mesh, E::Float64, nu::Float64, α::Float64)
    n    = length(mesh.points)
    ndim = mesh.ndim
    Fbc  = zeros(n*ndim)

    # Calculate average volumes
    #V    = [ cell_metric(c) for c in mesh.cells ]
    #Vavg = zeros(length(mesh.cells))
    #neig = get_neighbors(mesh.cells) # get list of neighbors per cell
    #for (i, cells) in enumerate(neig)
        #vol = V[i]
        #for c in cells
            #vol += V[c.id]
        #end
        #Vavg[i] = vol/(length(cells)+1)
    #end

    # volume average
    #V = 0.0
    #for c in mesh.cells
        #V += cellvolume(c)
    #end
    #V /= length(mesh.cells)

    for c in mesh.cells
        # get coordinates matrix
        np = length(c.points)
        C0 = cellcoords(c)

        #@show C0

        V  = cell_metric(c)
        #V = Vavg[c.id]
        #@show V

        T = cell_orientation(c)
        factor = V^(1./ndim)
        BC = basic_coords(c.shape)*factor

        #@show BC

        C = BC*T'

        # align C with cell orientation
        R, d = rigid_transform(C, C0)
        D = repmat( d , np, 1)
        C1 = C*R' + D
        U  = C1 - C0 # displacements matrix
        #@show U

        U  = vec(U')  # displacements vector

        K = matrixK(c, ndim, E, nu)
        #F = K*U*(1.-c.quality)
        #F = K*U*mesh.qmin/c.quality
        β = 0.97
        if mesh.qmin>β
            F  = K*U
        else
            F = K*U*(β-c.quality)/(β-mesh.qmin)*α
        end
        #F = K*U*(1.-c.quality^0.5)/(1.-mesh.qmin^0.5)*α
        #F  = K*U
        #exit()
        #F  = K*U*(1.-c.quality)/(1.-mesh.qmin)^0.2

        #FF = reshape(F, ndim, div(length(F),ndim))'

        # add forces to Fbc
        for (i,point) in enumerate(c.points)
            pid = point.id
            for j = 1:ndim
                Fbc[(pid-1)*ndim+j] += F[(i-1)*ndim+j]
            end
        end
    end

    return Fbc
end

function smooth!(mesh::Mesh; verbose=true, alpha::Float64=1.0, fixed::Bool=false, maxit::Int64=100, epsmin::Float64=1e-3,
    eps::Float64=1e-4, save_steps::Bool=false, filekey::String="smooth")

    verbose && println(BOLD, CYAN, "Mesh smoothing:", DEFAULT)

    E  = 1.0
    nu = 0.0

    ndim = mesh.ndim
    save_steps && save(mesh, "$filekey-0.vtk", verbose=false)

    Q = Float64[ c.quality for c in mesh.cells]
    q    = mesh.quality
    qmin = minimum(Q)
    dev  = stdm(Q, q)
    verbose && println("  histogram: ", hist(Q, 0.5:0.05:1.0)[2])
    verbose && @printf("  it: %3d  qmin: %7.5f  qavg: %7.5f  dev: %7.5f\n", 0, qmin, q, dev)

    for i=1:maxit

        # Lagrange multipliers matrix
        A   = mountA(mesh, fixed) 
        nbc = size(A,1)

        # Forces vector needed for smoothing
        F   = vcat( force_bc(mesh, E, nu, alpha), zeros(nbc) )

        # global stiffness plus LM
        K = mountKg(mesh, E, nu, A)

        # Solve system
        #@time U = K\F
        LUf = lufact(K)
        U = LUf\F

        # Update mesh
        for p in mesh.points
            p.x += U[(p.id-1)*ndim+1]
            p.y += U[(p.id-1)*ndim+2]
            if ndim>2
                p.z += U[(p.id-1)*ndim+3]
            end
        end

        for c in mesh.cells
            update!(c)
        end

        Q = Float64[ c.quality for c in mesh.cells]
        mesh.quality = sum(Q)/length(mesh.cells)
        mesh.qmin    = minimum(Q)
        save_steps && save(mesh, "$filekey-$i.vtk", verbose=false)

        if any(Q .== 0.0)
            error("smooth!: Smooth procedure got invalid element(s). Try using alpha<1.")
        end

        temp  = minimum(Q)
        Δq    = abs(q - mesh.quality)
        Δqmin = abs(qmin - temp)

        q    = mesh.quality
        qmin = temp
        dev  = stdm(Q, q)

        if Δq<eps && Δqmin<epsmin
            break
        end
        verbose && @printf("  it: %3d  qmin: %7.5f  qavg: %7.5f  dev: %7.5f", i, qmin, q, dev)
        verbose && println("  histogram: ", hist(Q, 0.5:0.05:1.0)[2])
    end

    verbose && println("  done.")

    return nothing
end

precompile(smooth!, (Mesh,) )

function laplacian_smooth(mesh::Mesh)
    # find nodal shares

    for i=1:length(mesh.cells)
        p = sum(..)/n
    end
end
