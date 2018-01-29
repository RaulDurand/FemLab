export solver2ex!

const Vect=Array{Float64,1}


mutable struct MecSolverData
    dom::Domain
    ips::Array{Ip,1}
    tol::Float64
    verbose::Bool
    umap::Array{Int,1} # list of unknown dofs
    pmap::Array{Int,1} # list of known dofs
    ndofs::Int         # number of dofs
end

# Assemble the global stiffness matrix (sparse matrix)
function mount_K(sdata::MecSolverData)

    R, C, V = Int64[], Int64[], Float64[]

    for elem in sdata.dom.elems
        Ke  = elem_stiffness(elem.mat, elem)
        map = elem_map(elem.mat, elem)
        nr, nc = size(Ke)
        for i=1:nr
            for j=1:nc
                push!(R, map[i])
                push!(C, map[j])
                push!(V, Ke[i,j])
            end
        end
    end

    return sparse(R, C, V, sdata.ndofs, sdata.ndofs)
end


function mount_RHS(dom::Domain, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in dom.elems
        map      = elem_map(elem.mat, elem)
        RHS[map] = elem_RHS(elem.mat, elem)
    end
    return RHS
end


function solve_inc(sdata::MecSolverData, K::SparseMatrixCSC{Float64, Int}, DU::Vect, DF::Vect)
    #  [  K11   K12 ]  [ U1? ]    [ F1  ]
    #  |            |  |     | =  |     |
    #  [  K21   K22 ]  [ U2  ]    [ F2? ]

    umap = sdata.umap
    pmap = sdata.pmap

    # Global stifness matrix
    ndofs = length(DU)
    nu    = length(umap)
    if nu == ndofs 
        warn("solve!: No essential boundary conditions.")
    end

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
    F2 = K22*U2
    U1 = zeros(nu)
    if nu>0
        RHS = F1 - K12*U2
        try
            LUfact = lufact(K11)
            U1  = LUfact\RHS
            F2 += K21*U1
        catch err
            warn("solve!: $err")
            U1 .= NaN
        end
    end

    # Completing vectors
    DU[1:nu]     = U1
    DF[nu+1:end] = F2
end

function solve_update_step(sdata::MecSolverData, K, ΔUa::Vect, ΔUi::Vect, R::Vect)
    ndofs = sdata.ndofs
    elems = sdata.dom.elems
    ips   = sdata.ips

    sdata.verbose && print("    solving...   \r")
    solve_inc(sdata, K, ΔUi, R)   # Changes unknown positions in ΔUi and R

    sdata.verbose && print("    updating... \r")

    # Restore the state to last converged increment
    for ip in ips; ip.data = deepcopy(ip.data0) end

    # Get internal forces and update data at integration points (update ΔFin)
    ΔFin = zeros(ndofs)
    ΔUt  = ΔUa + ΔUi
    for elem in elems  
        map = elem_map(elem.mat, elem)
        dU  = ΔUt[map]
        dF = elem_dF!(elem.mat, elem, dU) # gets internal force
        ΔFin[map] += dF
    end

    return ΔFin, R
end


function small_solver!(dom::Domain; nincs=1, maxits::Int=5, auto::Bool=false, NR::Bool=true, 
    tol::Number=1e-2, reset_bc::Bool=true, verbose::Bool=true, autosave::Bool=false, nout::Int=0)::Bool
    
    !NR      && (nincs = 1)
    (autosave && nout==0) && (nout=10)  # default value for nout
    autosave = nout>0

    if verbose
        print_with_color(:cyan,"FEM analysis:\n") 
        tic()
    end

    # get list of active elements
    aelems = [ elem for elem in dom.elems if elem.active ]

    # check if all active elements have material defined
    count = 0
    for elem in aelems
        if !isdefined(elem, :mat)
            count += 1
        end
    end
    count > 0 && error("There are $count active elements without material definition\n")

    # Calculate and apply boundary conditions
    for bc in dom.bcs
        apply_bc(bc)
    end

    # Arrays of prescribed dofs and unknown dofs
    dofs  = [dof for node in dom.nodes for dof in node.dofs]
    presc = [dof.prescU for dof in dofs]
    pdofs = dofs[presc]  # known dofs
    udofs = dofs[.!presc] # unknown dofs

    # Redefine tha array of dofs to list first the unknown ones
    dofs  = vcat(udofs, pdofs)
    ndofs = length(dofs)

    # Set equation id for each dof
    for (i,dof) in enumerate(dofs)
        dof.eq_id = i
    end
    
    # Fill arrays of prescribed dofs and unknown dofs
    umap  = [dof.eq_id for dof in udofs]
    pmap  = [dof.eq_id for dof in pdofs]

    # Get array with all integration points
    ips = [ ip for elem in aelems for ip in elem.ips ]

    # Global RHS vector 
    RHS = mount_RHS(dom, ndofs, 0.0)

#################bRYAN eXPRECOES###################################################################

    
    #Function definition for de natural and essential boundary conditions for alls dofs
    bryFun = Array{Any,1}(zeros(ndofs))
    
    for dof in dofs
        a = dof.bryF
        b = Expr(:(->),:(t,),a)
        bryFun[dof.eq_id] = eval(b)
    end

    F=Array{Float64,1}(zeros(ndofs))

    for i=1:ndofs
        F[i] = Base.invokelatest(bryFun[i], 0)
    end

######################################################################################################


    # Global U and F vectors
    U  = [ dof.bryU for dof in dofs]
    #F  = [ dof.bryF for dof in dofs] # nodal and face boundary conditions
    F += RHS

    # Solving process
    nu  = length(udofs)  # number of unknowns
    if verbose; println("  unknown dofs: $nu") end

    update_monitors!(dom)    # Tracking nodes, ips, elements, etc.

    if dom.nincs == 0 && autosave 
        save(dom, dom.filekey * "-0.vtk", verbose=false) # saves current domain state
        verbose && print_with_color(:green, "  $(dom.filekey)-0.vtk file written (Domain)\n")
    end

    # Set the last converged state at ips
    for ip in ips
        ip.data0 = deepcopy(ip.data)
    end

    # Incremental analysis
    T   = 0.0         # Pseudo time
    dT  = 1.0/nincs   # time increment

    dTs = 1.0/nout   # increment for saving vtk file
    Ts  = dTs        # T for next vtk file saving

    μdT = 1e-9       # minimum dT
    inc = 1          # current increment
    iout  = dom.nouts  # number of output file
    Fin   = zeros(ndofs)
    maxF  = norm(F) # Maximum norm of vectors of internal/external forces
    sdata = MecSolverData(dom, ips, tol, verbose, umap, pmap, ndofs) # solver data structure

    while T < 1.0 - μdT
        if verbose; print_with_color(:blue, "  increment $inc from T=$(round(T,10)) to T=$(round(T+dT,10)) (dT=$(round(dT,10))):\n") end
        ΔU, ΔF = dT*U, dT*F     # increment vectors
        R      = copy(ΔF)       # residual
        local ΔFin              # internal forces vector for current increment

        ΔUa    = zeros(ndofs)   # accumulated essential values (e.g. displacements)
        ΔUi    = copy(ΔU)       # values at iteration i
        nbigger  = 0            # counter to test non convergence

        # Newton Rapshon iterations
        residue   = 0.0
        converged = false
        maxfails  = 3   # maximum number of it. fails with close residues (> 0.9)
        nfails    = 0
        for it=1:maxits
            if it>1; ΔUi .= 0.0 end

            # Solve
            lastres = residue

            verbose && print("    assembling... \r")
            K = mount_K(sdata)

            verbose && print("    solving...   \r")
            solve_inc(sdata, K, ΔUi, R)   # Changes unknown positions in ΔUi and R

            verbose && print("    updating... \r")

            # Restore the state to last converged increment
            for ip in ips; ip.data = deepcopy(ip.data0) end

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin = zeros(ndofs)
            ΔUt  = ΔUa + ΔUi
            for elem in dom.elems  
                map = elem_map(elem.mat, elem)
                dU  = ΔUt[map]
                dF = elem_dF!(elem.mat, elem, dU) # gets internal force
                ΔFin[map] += dF
            end  # end of update

            maxF = max(norm(Fin+ΔFin), maxF)
            residue = maximum(abs, (ΔF-ΔFin)[umap] ) 

            # Update external forces vector
            # Update accumulated displacement
            ΔUa += ΔUi

            # Residual vector for next iteration
            R = ΔF - ΔFin  
            R[pmap] .= 0.0  # Zero at prescribed positions

            if verbose
                print_with_color(:bold, "    it $it  ")
                @printf("residue: %-10.4e", residue)
                println()
            end

            if residue < tol;        converged = true ; break end
            if isnan(residue);       converged = false; break end
            if it > maxits;          converged = false; break end
            if residue > 0.9*lastres;  nfails += 1 end
            if nfails == maxfails;     converged = false; break end
              
        end

        if converged

            Fin += ΔFin
            maxF = max(norm(Fin), maxF)

            # Store converged state at ips
            for ip in ips ip.data0 = deepcopy(ip.data) end

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.U += ΔUa[i]
                dof.F += ΔFin[i]
            end

            update_monitors!(dom) # Tracking nodes, ips, elements, etc.
            Tn = T + dT
            if Tn>=Ts && autosave
                iout += 1
                save(dom, dom.filekey * "-$iout.vtk", verbose=false, save_ips=save_ips)
                Ts = Tn - mod(Tn, dTs) + dTs
                verbose && print_with_color(:green, "  ", dom.filekey * "-$iout.vtk file written (Domain)\n")
            end

            inc += 1
            T   += dT
            if auto; dT = min(1.5*dT, 1.0/nincs, 1.0-T) end
        else
            if auto
                verbose && println("    increment failed.")
                dT *= 0.5
                if dT < μdT
                    print_with_color(:red, "solve!: solver did not converge\n",)
                    return false
                end
            else
                print_with_color(:red, "solve!: solver did not converge\n",)
                return false
            end
        end
    end

    # time spent
    if verbose
        spent = round(toq(),3)
        alls  = convert(Int, floor(Int, spent))
        msecs = convert(Int, round((spent-alls)*1000))
        hs    = div(alls, 3600)
        mins  = div(alls % 3600, 60)
        secs  = (alls % 3600) % 60
        println("  time spent: $(hs)h $(mins)m $secs.$(msecs)s")
    end
    
# Reset boudary conditions
    if reset_bc
        for node in dom.nodes
            for dof in node.dofs
                dof.bryU = 0.0
                dof.bryF = 0.0
                #dof.prescU = false
            end
        end
    end

    # Update number of used increments at domain
    dom.nincs += inc
    dom.nouts = iout

    return true

end

