using FemLab

# mesh 
bl = Block3D( [0 0 0; 1 1 0.5], nx=2, ny=2, nz=2)
mesh = generate_mesh(bl, verbose=false)

# fem domain
dom = Domain(mesh)
set_mat(dom.elems, DruckerPrager(E=100., nu=0.25, alpha=0.05, kappa=0.1) )
ip_dat = track(dom, dom.elems[1].ips[1])

# boundary conditions
set_bc( dom.nodes[:(z==0.0)] , ux=0, uy=0, uz=0)
set_bc( dom.nodes[:(z==0.5)] , uz=-0.033)
set_bc( dom.nodes[:(x==0 || x==1.0)], ux=0, uy=0)
set_bc( dom.nodes[:(y==0 || y==1.0)], ux=0, uy=0)

solve!(dom, nincs=10, verbose=true)

# boundary conditions
set_bc( dom.nodes[:(z==0.0)] , ux=0, uy=0, uz=0)
set_bc( dom.nodes[:(z==0.5)] , uz=0.008)
set_bc( dom.nodes[:(x==0 || x==1.0)], ux=0, uy=0)
set_bc( dom.nodes[:(y==0 || y==1.0)], ux=0, uy=0)

solve!(dom, nincs=10, verbose=true)

