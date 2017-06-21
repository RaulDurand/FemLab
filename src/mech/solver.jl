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

mutable struct MecSolverData
    dom::Domain
    ips::Array{Ip,1}
    tol::Float64
    verbose::Bool
    umap::Array{Int,1}
    pmap::Array{Int,1}
    ndofs::Int
end


function mount_K(sdata::MecSolverData)
    sdata.verbose && print("    assembling... \r")

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

    local S
    try
        S = sparse(R, C, V, sdata.ndofs, sdata.ndofs)
    catch err
        @show ndofs
        @show err
    end

    return S
end


function mount_RHS(dom::Domain, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in dom.elems
        map      = elem_map(elem.mat, elem)
        #RHS[map] = elem_RHS(elem, Δt::Float64)
        RHS[map] = elem_RHS(elem.mat, elem)
    end
    return RHS
end


function calc_residue(ΔF, ΔFin, Fmax, umap)
    abstol = 1e-8
    num    = norm((ΔF-ΔFin)[umap])
    #@show num
    #den    = max(norm(Fin), norm(F))
    den = Fmax
    #@show den
    if num < abstol 
        residue = num
    else
        residue = num/den
    end
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
            #U1  = K11\RHS
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


function solve!(dom::Domain; nincs=1, maxits::Int=5, auto::Bool=false, NR::Bool=true, scheme::String="FE", precision::Number=0.01,
    tol::Number=0.0, reset_bc::Bool=true, verbose::Bool=true, autosave::Bool=false, saveincs::Bool=false, nout::Int=0,
    save_ips::Bool=false)::Bool
    
    autosave && (saveincs = true)
    !NR      && (nincs = 1)
    (saveincs && nout==0) && (nout=10)  # default value for nout
    saveincs = nout>0

    if verbose
        print_with_color(:cyan,"FEM analysis:\n", bold=true) 
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
    pdofs = dofs[presc]
    udofs = dofs[.!presc]

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

    # Global U and F vectors
    U  = [ dof.bryU for dof in dofs]
    F  = [ dof.bryF for dof in dofs] # nodal and face boundary conditions
    F += RHS

    # Solving process
    nu  = length(udofs)  # number of unknowns
    if verbose; println("  unknown dofs: $nu") end

    update_monitors!(dom)    # Tracking nodes, ips, elements, etc.

    if dom.nincs == 0 && saveincs 
        save(dom, dom.filekey * "-0.vtk", verbose=false, save_ips=save_ips)
        verbose && print_with_color(:green, "  $(dom.filekey)-0.vtk file written (Domain)\n")
    end

    # Set the last converged state at ips
    for ip in ips
        ip.data0 = deepcopy(ip.data)
    end

    # State back-up for iterations
    data_bk = [ deepcopy(ip.data) for ip in ips ] 

    # Incremental analysis
    T   = 0.0
    dT  = 1.0/nincs

    dTs = 1.0/nout   # increment for saving vtk file
    Ts  = dTs        # T for next vtk file saving

    μdT = 1e-9       # minimum dT
    inc = 1
    iout = dom.nouts
    tol == 0.0 && (tol = precision)
    Fin   =  zeros(ndofs)
    maxF  = norm(F) # Maximum norm of vectors of internal/external forces
    sdata = MecSolverData(dom, ips, tol, verbose, umap, pmap, ndofs)

    while T < 1.0 - μdT
        if verbose; print_with_color(:blue, "  increment $inc from T=$(round(T,5)) to T=$(round(T+dT,5)) (dT=$(round(dT,10))):\n", bold=true) end  # color 111
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

            # Try FE step
            K = mount_K(sdata)
            ΔFin, R = solve_update_step(sdata, K, ΔUa, ΔUi, R)
            maxF = max(norm(Fin+ΔFin), maxF)
            #residue = calc_residue(ΔF, ΔFin, maxF, umap)
            residue = maximum(abs, (ΔF-ΔFin)[umap] ) #*****

            # use ME scheme
            if residue > tol && scheme == "ME"
                ΔFin, R = solve_step_ME(sdata, K, ΔUa, ΔUi, R)
                maxF = max(norm(Fin+ΔFin), maxF)
                #residue = calc_residue(ΔF, ΔFin, maxF, umap)
                residue = maximum(abs, (ΔF-ΔFin)[umap] ) #*****
            end

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
            if Tn>=Ts && saveincs
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

    if reset_bc
        for node in dom.nodes
            for dof in node.dofs
                dof.bryU = 0.0
                dof.bryF = 0.0
                dof.prescU = false
            end
        end
    end

    # Update number of used increments at domain
    dom.nincs += inc
    dom.nouts = iout

    return true

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

function solve_step_ME(sdata::MecSolverData, K1, ΔUa::Vect, ΔUi::Vect, R::Vect)
    # sdata:: in/out
    # ΔUa:: in
    # ΔUi:: in/out
    # R:: in/out

    # ME Corrector step:
    K2 = mount_K(sdata)
    K  = 0.5(K1 + K2)

    ΔFin, R = solve_update_step(sdata, K, ΔUa, ΔUi, R)

    return ΔFin, R
end

