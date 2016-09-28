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

export solve!

function mount_K(dom::Domain, ndofs::Int64)
    R, C, V = Int64[], Int64[], Float64[]

    for elem in dom.elems
        Ke  = elem_jacobian(elem)
        map = get_map(elem)
        nr, nc = size(Ke)
        for i=1:nr
            for j=1:nc
                push!(R, map[i])
                push!(C, map[j])
                push!(V, Ke[i,j])
            end
        end
    end
    
    return sparse(R, C, V, ndofs, ndofs)
end

function mount_RHS(dom::Domain, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in dom.elems
        map      = get_map(elem)
        #RHS[map] = elem_RHS(elem, Δt::Float64)
        RHS[map] = elem_RHS(elem)
    end
    return RHS
end

function solve!(dom::Domain; nincs::Int=1, maxits::Int=50, scheme::AbstractString="FE", precision::Float64=0.01, reset_bc::Bool=true, verbose::Bool=true, autosave::Bool=false, save_ips::Bool=false)

    if verbose; pbcolor(:cyan,"FEM analysis:\n") end

    # check if all elements have material defined
    count = 0
    for elem in dom.elems
        if !isdefined(elem, :mat)
            count += 1
        end
    end
    count > 0 && error("There are $count elements without material definition\n")

    # Set boundary conditions
    for bc in dom.node_bcs
        set_bc(bc.nodes; bc.conds...)
    end

    for bc in dom.face_bcs
        set_bc(bc.faces; bc.conds...)
    end

    for bc in dom.edge_bcs
        set_bc(bc.edges; bc.conds...)
    end

    # Fill array of dofs
    udofs = Array(Dof, 0)
    pdofs = Array(Dof, 0)

    for node in dom.nodes
        for dof in node.dofs
            if dof.prescU
                push!(pdofs, dof) 
            else
                push!(udofs, dof) 
            end
        end
    end

    dofs  = vcat(udofs, pdofs)
    ndofs = length(dofs)

    # Set equation id
    for (i,dof) in enumerate(dofs)
        dof.eq_id = i
    end
    
    # Fill arrays of prescribed dofs and unknown dofs
    presc = [ dof.prescU for dof in dofs ]
    pdofs = dofs[ presc]
    udofs = dofs[!presc]
    umap  = [dof.eq_id for dof in udofs]
    pmap  = [dof.eq_id for dof in pdofs]

    # Global RHS vector 
    RHS = mount_RHS(dom::Domain, ndofs::Int64, 0.0)

    # Global U F vectors
    U  = [ dof.bryU for dof in dofs]
    F  = [ dof.bryF for dof in dofs] # nodal and face boundary conditions
    F += RHS

    # Solving process
    nu  = length(udofs)
    if verbose; println("  unknowns dofs: $nu") end
    lam = 1.0/nincs

    DU, DF  = lam*U, lam*F
    residue = 0.0
    remountK = true
    tracking(dom) # Tracking nodes, ips, elements, etc.

    if autosave
        save(dom, dom.filekey * "-0" * ".vtk", verbose=false, save_ips=save_ips)
    end

    for inc=1:nincs
        if verbose; pcolor(:blue, "  increment $inc/$nincs:\n") end
        DU, DF = lam*U, lam*F
        R      = copy(DF) # residual

        DFa    = zeros(ndofs)
        DUa    = zeros(ndofs)
        nbigger= 0
        remountK = true

        converged = false
        for it=1:maxits
            if it>1; DU = zeros(ndofs) end
            solve_inc(dom, DU, R, umap, pmap, remountK, verbose)   # Changes DU and R
            #if it>1; DU[pmap] = 0.0 end

            if verbose; print("    updating... \r") end
            DUa += DU
            DFin = update!(dom.elems, dofs, DU) # Internal forces (DU+DUaccum?)

            R    = R - DFin
            DFa += DFin
        
            #residue = norm(R)
            lastres = residue
            residue = maxabs(R)

            if verbose
                pcolor(:bold, "    it $it  ")
                @printf("residue: %-15.4e", residue)
                println()
            end

            # Warning: use of remountK is bug prone!
            #if residue < 1e-10; remountK = false end
            #if residue < 1e-10 && it==1; remountK = false end
            #if it>1; remountK=true end

            if residue > lastres; nbigger+=1 end
            if residue<precision; converged = true ; break end
            if nbigger>15;        converged = false; break end
            if isnan(residue);    converged = false; break end
        end

        if converged
            tracking(dom) # Tracking nodes, ips, elements, etc.
            autosave && save(dom, dom.filekey * "-$inc" * ".vtk", verbose=false, save_ips=save_ips)
        end

        if !converged
            pcolor(:red, "solve!: solver did not converge\n",)
            return false
        end
    end

    if verbose && autosave
        pcolor(:green, "  $(dom.filekey)..vtk files written (Domain)\n")
    end

    if reset_bc
        for node in dom.nodes
            for dof in node.dofs
                dof.bryU = 0.0
                dof.bryF = 0.0
                dof.prescU = false
            end
        end
    end

    return true

end

function solve2!(dom::Domain; nincs::Int=1, maxits::Int=50, scheme::AbstractString="FE", precision::Float64=0.01, reset_bc::Bool=true, verbose::Bool=true, autosave::Bool=false, save_ips::Bool=false)
    
    if verbose; pbcolor(:cyan,"FEM analysis:\n") end

    # check if all elements have material defined
    count = 0
    for elem in dom.elems
        if !isdefined(elem, :mat)
            count += 1
        end
    end
    count > 0 && error("There are $count elements without material definition\n")

    # Set boundary conditions
    for bc in dom.node_bcs
        set_bc(bc.nodes; bc.conds...)
    end

    for bc in dom.face_bcs
        set_bc(bc.faces; bc.conds...)
    end

    for bc in dom.edge_bcs
        set_bc(bc.edges; bc.conds...)
    end

    # Fill array of dofs
    udofs = Array(Dof, 0)
    pdofs = Array(Dof, 0)

    for node in dom.nodes
        for dof in node.dofs
            if dof.prescU
                push!(pdofs, dof) 
            else
                push!(udofs, dof) 
            end
        end
    end

    dofs  = vcat(udofs, pdofs)
    ndofs = length(dofs)

    # Set equation id
    for (i,dof) in enumerate(dofs)
        dof.eq_id = i
    end
    
    # Fill arrays of prescribed dofs and unknown dofs
    presc = [ dof.prescU for dof in dofs ]
    pdofs = dofs[ presc]
    udofs = dofs[!presc]
    umap  = [dof.eq_id for dof in udofs]
    pmap  = [dof.eq_id for dof in pdofs]

    # Global RHS vector 
    RHS = mount_RHS(dom::Domain, ndofs::Int64, 0.0)

    # Global U F vectors
    U  = [ dof.bryU for dof in dofs]
    F  = [ dof.bryF for dof in dofs] # nodal and face boundary conditions
    F += RHS

    # Solving process
    nu  = length(udofs)
    if verbose; println("  unknowns dofs: $nu") end
    lam = 1.0/nincs

    DU, DF  = lam*U, lam*F
    residue = 0.0
    remountK = true
    tracking(dom) # Tracking nodes, ips, elements, etc.

    if autosave
        save(dom, dom.filekey * "-0" * ".vtk", verbose=false, save_ips=save_ips)
    end

    for inc=1:nincs
        if verbose; pcolor(:blue, "  increment $inc/$nincs:\n") end
        DU, DF = lam*U, lam*F
        R      = copy(DF) # residual

        #DFa    = zeros(ndofs)
        DUa    = zeros(ndofs)
        nbigger= 0
        remountK = true

        converged = false
        for it=1:maxits
            if it>1; DU = zeros(ndofs) end
            solve_inc(dom, DU, R, umap, pmap, remountK, verbose)   # Changes DU and R

            if verbose; print("    updating... \r") end
            DUa += DU
            DFin = update!(dom.elems, dofs, DUa) # Internal forces (DU+DUaccum?)

            DF[pmap] = DFin[pmap]  # Updates reactions
            R    = DF - DFin
            #DFa += DFin
        
            #residue = norm(R)
            lastres = residue
            residue = maxabs(R)

            if verbose
                pcolor(:bold, "    it $it  ")
                @printf("residue: %-15.4e", residue)
                println()
            end

            # Warning: use of remountK is bug prone!
            #if residue < 1e-10; remountK = false end
            #if residue < 1e-10 && it==1; remountK = false end
            #if it>1; remountK=true end

            if residue > lastres; nbigger+=1 end
            if residue<precision; converged = true ; break end
            if nbigger>15;        converged = false; break end
            if isnan(residue);    converged = false; break end
        end

        if converged
            tracking(dom) # Tracking nodes, ips, elements, etc.
            autosave && save(dom, dom.filekey * "-$inc" * ".vtk", verbose=false, save_ips=save_ips)
        end

        if !converged
            pcolor(:red, "solve!: solver did not converge\n",)
            return false
        end
    end

    if verbose && autosave
        pcolor(:green, "  $(dom.filekey)..vtk files written (Domain)\n")
    end

    # update state at ips
    for elem in dom.elems
        for ip in elem.ips
            ip.data0 = deepcopy(ip.data)
        end
    end

    if reset_bc
        for node in dom.nodes
            for dof in node.dofs
                dof.bryU = 0.0
                dof.bryF = 0.0
                dof.prescU = false
            end
        end
    end

    return true

end

# function solve_inc with static variables
let
    global solve_inc
    local  LUfact, K, K11, K12, K21, K22

function solve_inc(dom::Domain, DU::Vect, DF::Vect, umap::Array{Int,1}, pmap::Array{Int,1}, remountK::Bool, verbose::Bool)
    #  [  K11   K12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  K21   K22 ]  [ U2  ]    [ F2? ]

    # Global stifness matrix
    ndofs = length(DU)
    nu = length(umap)
    if nu == ndofs 
        pcolor(:red, "solve!: Warnig, no essential boundary conditions.\n")
    end

    if verbose; print("    assembling... \r") end

    if remountK
        K = mount_K(dom, ndofs)
        if nu>0
            nu1 = nu+1
            K11 = K[1:nu, 1:nu]
            K12 = K[1:nu, nu1:end]
            K21 = K[nu1:end, 1:nu]
        end
        K22 = K[nu+1:end, nu+1:end]
    end

    F1  = DF[1:nu]
    U2  = DU[nu+1:end]

    # Solve linear system
    if verbose; print("    solving...   \r") end
    F2 = K22*U2
    U1 = zeros(0)
    if nu>0
        RHS = F1 - K12*U2
        if remountK
            LUfact = lufact(K11)
            #@show cond(full(K11))
        end
        #U1  = K11\RHS
        U1 = LUfact\RHS
        F2 += K21*U1
    end

    # Completing vectors
    DU[1:nu]     = U1
    DF[nu+1:end] = F2
end

end # let


function update!(elems::Array{Element,1}, dofs::Array{Dof,1}, DU::Array{Float64,1})
    ndofs = length(dofs)
    DFin  = zeros(ndofs)

    # Update elements
    for elem in elems
        update!(elem, DU, DFin) # updates DFin
    end

    # Update dofs
    for (i,dof) in enumerate(dofs)
        dof.U += DU[i]
        dof.F += DFin[i]
    end

    return DFin

end

if VERSION >= v"0.4.0-dev+6521" 
    precompile(solve!,(Domain,))
end
