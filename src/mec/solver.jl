export solve!

function mount_K(dom::Domain, ndofs::Int64)
    R, C, V = Int64[], Int64[], Float64[]

    for elem in dom.elems
        Ke  = stiff(elem)
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
    
    sparse(R, C, V, ndofs, ndofs)
end

function solve!(dom::Domain; nincs::Int=1, scheme::String="FE", precision::Float64=0.01, reset_bc::Bool=true, verbose::Bool=true)
    if verbose; println("FEM analysis:") end
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

    dofs  = [udofs, pdofs]
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

    # Global U F vectors
    U = [ dof.bryU for dof in dofs]
    F = [ dof.bryF for dof in dofs]
    #println(U)

    # Solving process
    nu  = length(udofs)
    if verbose; println("  unknowns dofs: $nu") end
    lam = 1.0/nincs
    #println(lam)
    DU, DF  = lam*U, lam*F
    residue = 0.0

    for inc=1:nincs
        if verbose; println("  increment $inc:") end
        DU, DF = lam*U, lam*F
        R      = copy(DF) # residual
        #println(R)
        DFa    = zeros(ndofs)

        maxits    = 100
        converged = false
        for it=1:maxits
            if it>1; DU = zeros(ndofs) end
            solve_inc(dom, DU, R, umap, pmap, verbose)   # Changes DU and R
            if verbose; print("    updating... \r") end
            DFin = update!(dom.elems, dofs, DU) # Internal forces (DU+DUaccum?)

            R    = R - DFin
            DFa += DFin
        
            #residue = maxabs(DF[umap] - DFa[umap])
            residue = norm(R)
            tracking(dom) # Tracking nodes, ips, elements, etc.

            if verbose; println("    it $it  residue: $residue") end
            if residue<precision; converged = true ; break end
            if residue>100.;      converged = false; break end
        end

        if !converged
            error("solve!: solver did not converge")
        else
            continue
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

end


function solve_inc(dom::Domain, DU::Vect, DF::Vect, umap::Array{Int,1}, pmap::Array{Int,1}, verbose::Bool)
    #  [  K11   K12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  K21   K22 ]  [ U2  ]    [ F2? ]

    # Global stifness matrix
    ndofs = length(DU)
    nu = length(umap)
    if verbose; print("    assembling... \r") end
    K = mount_K(dom, ndofs)
    #println(full(K))
    if nu>0
        nu1 = nu+1
        K11 = K[1:nu, 1:nu]
        K12 = K[1:nu, nu1:end]
        K21 = K[nu1:end, 1:nu]
    end
    K22 = K[nu+1:end, nu+1:end]

    F1  = DF[1:nu]
    U2  = DU[nu+1:end]

    # Solve linear system
    if verbose; print("    solving...   \r") end
    F2 = K22*U2
    U1 = zeros(0)
    if nu>0
        RHS = F1 - K12*U2
        U1  = K11\RHS
        F2 += K21*U1
    end

    # Completing vectors
    #umap and pmap are not needed...
    DU[umap] = U1
    DF[pmap] = F2
end


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

