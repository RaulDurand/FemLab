using FemLab

coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]

blt = BlockTruss(coord, conn)
mesh = generate_mesh(blt, verbose=false)

dom = Domain(mesh)

set_mat(dom.elems, Truss(E=6.894757e7, A=0.043) )

set_bc( dom.nodes[:(x==0 && y==0)] , ux=0, uy=0)
set_bc( dom.nodes[:(x==0 && y==9)] , ux=0, uy=0)
set_bc( dom.nodes[:(x==9 && y==0)] , fy=-450.)
set_bc( dom.nodes[:(x==18&& y==0)] , fy=-450.)

truss = dom


function get_chrom(truss::Domain)
    # Get the chromosome (array) from an individual (truss) 
    return Float64[ elem.mat.A for elem in truss.elems ]
end

function mfactor()
    return 0.5 + rand()
end

function mutate(cr, rate=0.5)
    n = length(cr)
    # ngens: number of gens to mutate
    ngens = ifloor(n*rate)  # int or ifloor
    idx = [1:n]
    idx = shuffle(idx)[1:ngens]
    for i in idx
        cr[i] *= mfactor()
    end
    return cr
end

function crossover(cr1, cr2)
    # Generates a new chromosome based on two parents chromosomes
    n   = length(cr1)
    pos = rand(1:n-1)
    cr  = [ cr1[1:pos], cr2[pos+1:end] ]
    return cr
end


function gen_ind()
    # Generates a chromosome for a new individual
    cr = get_chrom(truss)
    return mutate(cr)
end


function gen_pop(num::Int)
    # Generates the population
    # gen_ind: function that generates an individual
    # num    : number of individuals in the population
    pop = [ gen_ind() for i=1:num ]
    return pop
end

function fitness(cr)
    # Calculates the fitness for a single individual (cromossome)

    n = length(cr) # equal to number of bars
    # reassign chromosome values to individuals
    E = truss.elems[1].mat.E
    for i=1:n
        #truss.elems[i].mat.A = cr[i]
        set_mat( truss.elems[i], Truss(E=E, A=cr[i] ))
    end

    # calculate weight
    weight = 0.0
    gamma  = 25.9  # unit weigth kN/m3
    for elem in truss.elems
        c1 = elem.nodes[1].X
        c2 = elem.nodes[2].X
        L  = norm(c2-c1)
        A  = elem.mat.A
        weight = weight + A*L*gamma
    end
    
    # reset stresses and strains
    reset(truss.nodes) # reset displacements
    reset(truss.elems) # reset stresses and strains

    # solving
    solve!(dom, verbose=false, reset_bc=false)

    # checking displacements
    #disps = [ abs(node[:uy].U) for node in truss.nodes ]
    disps = [ abs(node.dofdict[:uy].U) for node in truss.nodes ]
    max_disp = maximum(disps)
    if max_disp > 0.0508 # displacement (m)
        return -1. # flag
    end

    # checking stresses
    stresses = [ elem.ips[1].data.σ for elem in truss.elems] 
    tens_stress = abs(maximum(stresses))
    comp_stress = minimum(stresses)

    if tens_stress>130000.0 # stress (KPa)
        return -1
    end

    return 1./weight
end

function evolve(n::Int; nger::Int=5, elit::Float64=0.1, mutation::Float64=0.3)
    # Performs the evolution
    pop  = gen_pop(n)
    data = DTable( [ :generation, :fitness ] )

    for k=1:nger
        println("\nGeneration: ", k)

        fits = [ fitness(ind) for ind in pop ]
        idxs = sortperm(fits, rev=true)
        pop  = pop[idxs]      # orders pop according to fitness
        fits = fits[idxs]     # orders fits

        bestw = 1 / fits[1]
        println("Best fitness :", fits[1])
        println("Best Weight  :", bestw)

        push!(data, [k, fits[1] ] )
        
        if k==nger; break end

        nelit  = ifloor(elit*n)
        ncross = n-nelit
        pop = pop[1:nelit]

        # Cross over
        for i=1:ncross
            ind1 = pop[rand(1:nelit)]
            ind2 = pop[rand(1:nelit)]
            ind  = crossover(ind1, ind2)
            push!(pop, ind)
        end

        # Mutation
        nmutate = ifloor(mutation*n)
        for i=1:nmutate
            pos = rand(nelit+1:n)
            ind = pop[pos]
            mutate(ind)
        end

    end
    #print "\nFinal population:", pop
    save(data, "fitness.dat")
end

srand(0)
evolve(90, nger=100, elit=0.1, mutation=0.3)
#@time evolve(90, nger=10, elit=0.1, mutation=0.3)

readline()
