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

#Solids stiffness matrix for computing frequencies

function mount_Kf(sdata::MecSolverData)

    R, C, V = Int64[], Int64[], Float64[]
    elems=sdata.dom.elems[:solids]

    for elem in elems
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

# Assemble the global mass matrix (sparse matrix)

function mount_M(sdata::MecSolverData)

    R, C, V = Int64[], Int64[], Float64[]

    for elem in sdata.dom.elems
        Me  = elem_mass(elem.mat, elem)
        map = elem_map(elem.mat, elem)
        nr, nc = size(Me)
        for i=1:nr
            for j=1:nc
                push!(R, map[i])
                push!(C, map[j])
                push!(V, Me[i,j])
            end
        end
    end

    return sparse(R, C, V, sdata.ndofs, sdata.ndofs)
end

#Solids mass matrix for computing frequencies

function mount_Mf(sdata::MecSolverData)

    R, C, V = Int64[], Int64[], Float64[]
    elems=sdata.dom.elems[:solids]
 
    for elem in elems
        Me  = elem_mass(elem.mat, elem)
        map = elem_map(elem.mat, elem)
        nr, nc = size(Me)
        for i=1:nr
            for j=1:nc
                push!(R, map[i])
                push!(C, map[j])
                push!(V, Me[i,j])
            end
        end
    end
    
    return sparse(R, C, V, sdata.ndofs, sdata.ndofs)
end

# Two first frequencies for Rayleigh Ritz method of damping
function Frequency(sdata::MecSolverData, K::SparseMatrixCSC{Float64, Int}, M::SparseMatrixCSC{Float64, Int}, nmod::Int64)

    #  [  K11   K12 ]  [ U1? ] 
    #  |            |  |     | 
    #  [  K21   K22 ]  [ U2  ] 

    #  [  M11   M12 ]  [ U1? ] 
    #  |            |  |     | 
    #  [  M21   M22 ]  [ U2  ] 


    umap = sdata.umap #vetor com unknow
    pmap = sdata.pmap #vetor com know

    #Eigenvals problem solution K22 - w²*M22 = 0 ==> (M22^-1)*K22 - w²*I = 0 ==> P - w²*I = 0

    ndofs = size(K,1) 
    nu    = length(umap) 

    if nu == ndofs 
        warn("solve!: No essential boundary conditions.")
    end

    #Lumped Mass Matrix

    for i=1:size(M,1)
        for j=1:size(M,2)
            if i != j
                M[i,i] += M[i,j]
                M[i,j] = 0.0
            end
        end
    end


    if nu>0
        K11 = K[1:nu, 1:nu]
        M11 = M[1:nu, 1:nu]
    end

    #M11=full(M11)
    #K11=full(K11)

    #Delete columns and rows that corresponding at bars and joint elements

    nodesbars=get_nodes(sdata.dom.elems[:lines]) #nodes of bars

    idsd = Array{Int64,1}(0) #idsd ids for delete

    for node in nodesbars
        for dof in node.dofs
            push!(idsd,dof.eq_id)
        end
    end
    idsd = unique(idsd)
    sort!(idsd)
    nidsd = length(idsd)

    for i=1:nidsd
       M11=M11[1:size(M11,1)     .!= idsd[nidsd+1-i,1],:];
       M11=M11[:,1:size(M11,2) .!= idsd[nidsd+1-i,1]];
       K11=K11[1:size(K11,1)     .!= idsd[nidsd+1-i,1],:];
       K11=K11[:,1:size(K11,2) .!= idsd[nidsd+1-i,1]];
    end 

    #M11V = zeros(size(M11,2))  #Reduced Matrix Mass on vector form
    M11IV = zeros(size(M11,2)) #Vetor of inverse matrix mass 

    for i=1:size(M11,2)
        #M11V[i] = M11[i,i]
        M11IV[i] = 1/M11[i,i]
    end

    #M11IV = 1./diag(M11)

    #MI = inv(M11)
    #P = MI * K11
    P = M11IV.*K11
    #P = full(P)

    #Aulovals solution problem for P matrix

	#Eig=eigfact(P);
    nevn=size(P,2)-2 #eigvals number
    Eig=eigs(P, nev=20, which=:SM)


    #w0 = Eig[:values];
    w0 = Eig[1]
    wi = copy(w0)
    #Select Logic and possible vals

    wi=real(wi[isreal.(wi)]) #delete imaginarian numbers
    #wi=unique(wi) #all vals diferents
    wi=wi[wi[1:size(wi)[1]] .>0] # delete zero and negative vals
    wi=sort!(wi) # sort the vector
    
    w = Array{Float64,1}(0)

    for i=1:nmod
        push!(w,wi[i])
    end

   
    v = Array{Float64,2}(size(P,2),nmod)

    for i=1:nmod
        col = findfirst(w0,w[i])
        v[:,i] = Eig[2][:,col]
        #v[:,i] = Eig[:vectors][:,col]
    end
     
    w = w.^0.5
    @show w

	return w, v

end

#Rayleigh-Ritz method for damping

function DampingRayleigh(M,K,w,xi1,xi2)


    alfaD = 2 * ((w[1] ^ 2) * w[2] * xi2 - w[1] * (w[2] ^ 2) * xi1) / ((w[1] ^ 2)-(w[2] ^ 2));

    betaD= 2 * (w[1] * xi1 - w[2] * xi2)/((w[1] ^ 2) - (w[2] ^ 2));

    C = alfaD*M + betaD*K;

    return C

end

#Initial aceleration

function initialaceleration(sdata::MecSolverData, K::SparseMatrixCSC{Float64, Int}, M::SparseMatrixCSC{Float64, Int},
C::SparseMatrixCSC{Float64, Int}, U::Vect, V::Vect, A::Vect, F::Vect)

 #Dynamic equilibrium.(Ua0 = (M22^-1)*(F0 - C22Uv0 - K22Ud0)

    #  [  K11   K12 ]  [ U1  ]  [  C11   C12 ]  [ V1  ]   [  M11   M12 ]  [ A1? ]   [ F1  ]
    #  |            |  |     | +|            |  |     | + |            |  |     | = |     |
    #  [  K21   K22 ]  [ U2  ]  [  C21   C22 ]  [ V2  ]   [  M21   M22 ]  [ A2  ]   [ F2? ]

    umap = sdata.umap
    pmap = sdata.pmap

    ndofs = length(U)  
    nu    = length(umap) 

    if nu == ndofs 
        warn("solve!: No essential boundary conditions.")
    end

   if nu>0
        nu1 = nu+1

        K11 = K[1:nu, 1:nu]
        K12 = K[1:nu, nu1:end]
        K21 = K[nu1:end, 1:nu]

        C11 = C[1:nu, 1:nu]
        C12 = C[1:nu, nu1:end]
        C21 = C[nu1:end, 1:nu]

        M11 = M[1:nu, 1:nu]
        M12 = M[1:nu, nu1:end]
        M21 = M[nu1:end, 1:nu]

    end

    K22 = K[nu+1:end, nu+1:end]
    C22 = C[nu+1:end, nu+1:end]
    M22 = M[nu+1:end, nu+1:end]

    F1   = F[1:nu]
    U1  = U[1:nu]
    V1  = V[1:nu]
    U2  = U[nu+1:end]
    V2  = V[nu+1:end]
    A2  = A[nu+1:end]

    # Solve linear system

    F2 = K22*U2 + C22*V2 + M22*A2  
    A1 = zeros(nu)

    if nu>0
        RHS = F1 - K11*U1 - K12*U2 - C11*V1 - C12*V2 - M12*A2
        
        try
            LUfact = lufact(M11)
            A1  = LUfact\RHS
            F2 += K21*U1 + C21*V1 + M21*A1
        catch err
            warn("solve!: $err")
            A1 .= NaN
        end
    end

    # Completing vectors
    A[1:nu]     = A1
    F[nu+1:end]   = F2

end

#Newmark method to find a increment of displasment DU due a load increment DF

function solve_incnewmark(sdata::MecSolverData, K::SparseMatrixCSC{Float64, Int}, M::SparseMatrixCSC{Float64, Int},
C::SparseMatrixCSC{Float64, Int}, DU::Vect, DF::Vect, Da::Vect, Vt::Vect, At::Vect, Deltat::Float64)

#Pseudo-stiffness matrix

Kp = K + (4/(Deltat^2))*M + (2/Deltat)*C

#Pseudo-forces vector

DFp = DF + M*(At + 4*Vt/Deltat + 4*(-Da)/(Deltat^2)) + C*(Vt + (2*(-Da)/Deltat))

# Solve the Kp * DU = DFp system

    #  [  Kp11   Kp12 ]  [ DU1? ]    [ DFp1  ]
    #  |              |  |      | =  |       |
    #  [  Kp21   Kp22 ]  [ DU2  ]    [ DFp2? ]

    umap = sdata.umap
    pmap = sdata.pmap
    ndofs = length(DU)
    nu    = length(umap)

    if nu == ndofs 
        warn("solve!: No essential boundary conditions.")
    end

    if nu>0
        nu1 = nu+1
        Kp11 = Kp[1:nu, 1:nu]
        Kp12 = Kp[1:nu, nu1:end]
        Kp21 = Kp[nu1:end, 1:nu]
    end

    Kp22 = Kp[nu+1:end, nu+1:end]
    DFp1  = DFp[1:nu]
    DU2  = DU[nu+1:end]

    # Solve linear system

    DFp2 = Kp22*DU2
    DU1 = zeros(nu)

    if nu>0
        RHS = DFp1 - Kp12*DU2
        try
            LUfact = lufact(Kp11)
            DU1  = LUfact\RHS
            DFp2 += Kp21*DU1
        catch err
            warn("solve!: $err")
            U1 .= NaN
        end
    end

    # Completing vectors

    DU[1:nu]     = DU1
    DFp[nu+1:end] = DFp2
    DF = DFp - ( M*(At + 4*Vt/Deltat + 4*(-Da)/(Deltat^2)) + C*(Vt + (2*(-Da)/Deltat)))
       
end

function SismicForce(sdata::MecSolverData, M::SparseMatrixCSC{Float64, Int}, F::Vect, AS::Array{Float64,2}, keysis::Symbol,
tfia::Float64, tds::Float64)


    ndat = length(AS) #quantity of aceleration data
    AS2 = zeros(ndat+1) #Sismic aceleration data how vector
    c = 0
    for i=1:size(AS,1)
        for j=1:size(AS,2)
            c += 1
            AS2[c+1] = AS[i,j]
        end
    end

    vts = zeros(ndat+1) #time vetor correspond to acelerations
    
    for i=1:ndat+1
        vts[i] = (i-1)*tds/ndat
    end
    
    FAS = hcat(vts,AS2) # Function of aceleration
    #Interpolation of the aceleration value

        inic = 0
        fin = 0

        for i=1:ndat
            if FAS[i,1]<=tfia<=FAS[i+1,1]
                inic = i
                fin = i+1
            end
            if inic!=0 & fin!=0; break end
        end
    
        m = (FAS[fin,2]-FAS[inic,2])/(FAS[fin,1]-FAS[inic,1])
        acel = FAS[inic,2] + m*(tfia - FAS[inic,1])
    #@show acel

    #Dof sismic aceleration    

    VAS  = zeros(sdata.ndofs) #Sismic aceleration vector accord dof

    nodes = sdata.dom.nodes
    
    for node in nodes
        dof = node.dofdict[keysis]
        VAS[dof.eq_id] += acel
    end    
    
    #Dof sismic force
    FS = M*VAS
    F += FS

    return F

end

 




#Dynamic solver with newtown raphson method for non liner

function small_solverdyn!(dom::Domain, ta::Float64, nint::Int64, ndom::Int64=0, nmod::Int64=10, tds::Float64=0.0,
tss::Float64=0.0;owse=false, sism=false, nincs=1, maxits::Int=5, auto::Bool=false, NR::Bool=true, 
    tol::Number=1e-2, reset_bc::Bool=true, verbose::Bool=true, autosave::Bool=false, nout::Int=0)::Bool


    #If the problem is sismic, read the sismic acelerations asking to user the file's name

    if sism==true
         print("What is the .dat file name of the sismic acelerations?")
         AS = readdlm("$(chomp(readline())).dat")
         AS= 9.81*AS
         print("What is the key correspond to sismic direction (fx, fy, fz)?")
         keysis = parse(chomp(readline()))
    end        

    #ndom = dommains number. Quantity of outs domains that user want save excepting the end. <= nint.

    #Duration of interval time i Deltat = ta/nint

    Deltat = ta/nint

    tfia = 0.0 #tiempo final of actual interval


    !NR      && (nincs = 1)
    (autosave && nout==0) && (nout=10)  # default value for nout
    autosave = nout>0

    if verbose
        print_with_color(:cyan,"FEM analysis/:\n") 
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

    # MecSolverData

    sdata = MecSolverData(dom, ips, tol, verbose, umap, pmap, ndofs) # solver data structure #no sei como funciona

    #Global U, V, A and F vectors

    #Function definition for de natural and essential boundary conditions for alls dofs
    
    #Array function of natural bcs bryFun
    bryFun = Array{Any,1}(zeros(ndofs))
    
    for dof in dofs
        a = dof.bryF
        b = Expr(:(->),:(t,),a)
        bryFun[dof.eq_id] = eval(b)
    end

    Ui = [ dof.bryU for dof in dofs]  #Initial displacements Ui of the first interval
    Fi=Array{Float64,1}(zeros(ndofs)) #Initial forces of the first interval
    
    for i=1:ndofs
        Fi[i] = Base.invokelatest(bryFun[i], 0)
    end

    #If the problem has a sism, the force sismic is add

    if sism==true && tss<=tfia
        M = mount_M(sdata)
        Fi = SismicForce(sdata, M, Fi,AS,keysis,tfia,tds)
    end


    Ai=zeros(ndofs)
    
   #Initial displacement and velocity (nuls) 

    Ui=zeros(ndofs)
    Vi = zeros(ndofs)

#If static efects of own weight is accont (static efects of own weight owse owse =true). Make first static analise with
    #small_solver and weight loads.
    if owse == true
         Ui = [ dof.U for dof in dofs]
    end

   

    #M,K,C matrixs for initial acelerations
    
    #Mass matrix for frequencies computing (only solid elements)
    verbose && print("    assembling... \r")
    Mf = mount_Mf(sdata)
    #Stiffness matrix for frequencies computing
    verbose && print("    assembling... \r")
    Kf = mount_Kf(sdata)

    #Two first frequencies
           
    w1, v = Frequency(sdata, Kf, Mf, nmod)


    #Modals Deformed

    nodessolids = get_nodes(dom.elems[:solids]) #nodes of solids

    idss = Array{Int64,1}(0) #idss ids of dof solids

    for node in nodessolids
        for dof in node.dofs
            push!(idss,dof.eq_id)
        end
    end
    idss = unique(idss)
    sort!(idss)
    nidss = length(idss)
    idefmod = 1    

    for j=1:nmod

        Umod = zeros(ndofs)

        for i=1:length(v[:,j])
            Umod[idss[i]] = v[i,j]
        end

        for (k,dof) in enumerate(dofs)
            dof.U = Umod[k] 
        end

        save(dom, dom.filekey * "-defmod$idefmod.vtk", verbose=false) # saves current domain state for modal deformated
        verbose && print_with_color(:green, "  ", dom.filekey * "-defmod$idefmod.vtk file written (Domain)\n")
        idefmod += 1


    end
  
    verbose && print("    assembling... \r")
    M = mount_M(sdata)

    verbose && print("    assembling... \r")
    K = mount_K(sdata)
           
    #Rayleigh-Ritz method for Damping

    xi1=dom.elems[1].mat.xi
    xi2=dom.elems[1].mat.xi

    #Interacting with the user for select the correct frequêncies based on modal deformated

    w = zeros(2)
    println("See the vtk files of the modal deformated and decide")
    print("What is the first frequency that you want?")
    w[1] = w1[parse(chomp(readline()))]
    print("What is the second frequency that you want?")
    w[2] = w1[parse(chomp(readline()))]

    C=DampingRayleigh(M,K,w,xi1,xi2)

    #Initial acelerations

         #Initial Displacements and acelerations for dynamic analise
    
            U0=zeros(ndofs)
            V0=zeros(ndofs)

    initialaceleration(sdata,K,M,C,U0,V0,Ai,Fi)

    # Solving process
    nu  = length(udofs)  # number of unknowns
    if verbose; println("  unknown dofs: $nu") end

   

    for (i,dof) in enumerate(dofs)

        dof.U = Ui[i] 
        dof.V = Vi[i] 
        dof.A = Ai[i] 
        dof.F = Fi[i]
                                           
    end

    update_monitors!(dom)    # Tracking nodes, ips, elements, etc.

    if ndom != 0  
        save(dom, dom.filekey * "-0.vtk", verbose=false) # saves current domain state
        verbose && print_with_color(:green, "  $(dom.filekey)-0.vtk file written (Domain)\n")
    end

    #What interval of time save ? 4V solvedynlinear wint. Identification of saved domain idom.

    wint = round(nint/ndom)
    ctewint = copy(wint)
    idom = 1

    # Set the last converged state at ips
    for ip in ips
        ip.data0 = deepcopy(ip.data)
    end

    #Loop for iteration of time
    
    #Accumulated stiffnes force Fin
    Fin   = zeros(ndofs) # internal stiffness force
    lastFin = zeros(ndofs) # internal stiffness force of the last converged iteration

    for l=1:nint

        T   = 0.0         # Pseudo time
        dT  = 1.0   # time increment
        inc = 1 
        tiia = (l-1)*Deltat #initial time of actual interval

        while T < 1.0

            Fin = copy(lastFin)
            
            Deltata = dT*Deltat # actual time interval Deltata. This change with the automatic pass.
            tfia = tiia + Deltata #final time of actual interval

            #Final displacements and forces of the l interval  

            Uf = zeros(ndofs)
            Ff = zeros(ndofs)

            for i=1:ndofs
                Ff[i] = Base.invokelatest(bryFun[i], tfia)
            end

            #If the problem has a sism, the force sismic is add

            if sism==true && tss<=tfia && tds+tss>=tfia
                M = mount_M(sdata)
                Ff = SismicForce(sdata, M,Ff,AS,keysis,tfia,tds)
               
#= 
                if l==100
                    sss=0      
                    for node in dom.nodes
                        dofsss = node.dofdict[keysis]
                        sss += Ff[dofsss.eq_id] 
                    end     
                    @show sss
                    exit()
                end
=#            
            end

            if verbose; print_with_color(:blue, "tfia=$tfia s, increment $inc from T=$(round(T,10)) to T=$(round(T+dT,10)) (dT=$(round(dT,10))):\n") end
            ΔU, ΔF = (Uf), (Ff)     # increment vectors
            R      = copy(ΔF)       # residual
            local ΔFin              # internal forces vector for current increment

            ΔUa    = zeros(ndofs)   # accumulated essential values (e.g. displacements)
            ΔUi    = copy(ΔU)       # values at iteration i
            nbigger  = 0            # counter to test non convergence

            # Loop for Newton Rapshon iterations

            residue   = 0.0
            converged = false
            maxfails  = 3   # maximum number of it. fails with close residues (> 0.9)
            nfails    = 0

            Ait, Vit, TFin = zeros(ndofs), zeros(ndofs), zeros(ndofs)


            R -= Fin

            for it=1:maxits
                
                if it>1; ΔUi .= 0.0 end
    

                # Solve
                lastres = residue
                #Stiffness matrix
                verbose && print("    assembling... \r")
                K = mount_K(sdata)
                
                #Damping matrix
                C = DampingRayleigh(M,K,w,xi1,xi2)

                verbose && print("    solving...   \r")
                solve_incnewmark(sdata, K, M, C, ΔUi, R, ΔUa, Vi, Ai, Deltata)   # Changes unknown positions in ΔUi and R

                verbose && print("    updating... \r")
             
                # Restore the state to last converged increment
                for ip in ips; ip.data = deepcopy(ip.data0) end

                # Get internal forces and update data at integration points (update ΔFin)
                ΔFin = zeros(ndofs)
                ΔUt  = ΔUa + ΔUi #+ Ui #total displacement until now

                #@show ΔUt
                
                #Stiffness forces for total displacement
                for elem in dom.elems  
                    map = elem_map(elem.mat, elem)
                    dU  = ΔUt[map]
                    dF = elem_dF!(elem.mat, elem, dU) # gets internal force
                    ΔFin[map] += dF
                    #Fin[map] += dF
                end  # end of update

                Fin += ΔFin

                #Mass and damping forces adition for the total of internal forces of the present iteration
                Vit = -Vi + 2*(ΔUa + ΔUi)/Deltata;
                Ait = -Ai + 4*((ΔUa + ΔUi) - Vi*Deltata)/(Deltata^2);               
                #@show Ait

                #Total internal(stiffness), inertial and damping forces TFin 
                TFin = Fin + C*Vit + M*Ait

                residue = maximum(abs, (ΔF-TFin)[umap] )
                
                # Update external forces vector
                # Update accumulated displacement
                ΔUa += ΔUi

                # Residual vector for next iteration for the Newmark method
                R = ΔF - Fin  
                R[pmap] .= 0.0  # Zero at prescribed positions
                if verbose
                    print_with_color(:bold, "    it $it  ")
                    @printf("residue: %-10.4e", residue)
                    println()
                end

                #Fin -= ΔFin #get back to the initial internal force before of the NR iterations
                if residue > tol; Fin -= ΔFin end
                if residue < tol;        converged = true ; break end
                if isnan(residue);       converged = false; break end
                if it > maxits;          converged = false; break end
                if residue > 0.9*lastres;  nfails += 1 end
                if nfails == maxfails;     converged = false; break end
            end
            if converged
                
                lastFin = copy(Fin)

                # Store converged state at ips
                for ip in ips ip.data0 = deepcopy(ip.data) end

                # Update nodal variables at dofs of actual time interval

                Af = copy(Ait)
                Vf = copy(Vit)
                Uf = copy(ΔUa)+Ui
                Ff = copy(TFin)

                #Initial Displacements and Forces of the next interval l+1
    
                Ui=copy(Uf)
                Fi=copy(Ff)
                Ai=copy(Af)
                Vi=copy(Vf)

                if (T+dT)==1 #if corresponded of the final point of the non modifiqued interval 
                
                    for (i,dof) in enumerate(dofs)
                        dof.U = Uf[i] 
                        dof.V = Vf[i] 
                        dof.A = Af[i] 
                        dof.F = Ff[i]                                           
                    end

                    update_monitors!(dom) # Tracking nodes, ips, elements, etc.
                
                
                    if wint == l  
                        save(dom, dom.filekey * "-$idom.vtk", verbose=false) # saves current domain state
                        verbose && print_with_color(:green, "  ", dom.filekey * "-$idom.vtk file written (Domain)\n")
                        wint += ctewint
                        idom += 1
                    end
                end        
    
                T += dT
                tiia += Deltata 
                inc += 1

            else
                if auto
                    verbose && println("    increment failed.")
                    dT *= 0.5
                    if dT < 0.01*Deltat
                        print_with_color(:red, "solve!: solver did not converge\n",)
                        return false
                    end
                else
                    print_with_color(:red, "solve!: solver did not converge\n",)
                    return false
                end
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
                dof.bryF = :(0*t)
                dof.prescU = false
            end
        end
    end

    return true

end

