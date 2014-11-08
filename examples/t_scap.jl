using FemLab

# mesh 
bl = Block3D( [0 0 0; 1 1 0.5], nx=1, ny=1, nz=1)
mesh = generate_mesh(bl)

# fem domain
dom = Domain(mesh)
set_mat(dom.elems, SingleCap(E=100., nu=0.25, alpha=0.05, kappa=0.1, rho=8., H=0.) )
#set_mat(dom.elems, DruckerPrager(E=100., nu=0.25, alpha=0.05, kappa=0.1) )
ip_dat = track(dom, dom.elems[1].ips[1])

#dom.elems[1].ips[1].data.pl = true

# boundary conditions
set_bc( dom.nodes[:(z==0.0)] , ux=0, uy=0, uz=0)
#set_bc( dom.nodes[:(z==0.5)] , uz=-0.04)
set_bc( dom.faces[:(z==0.5)] , tz=-2.)
set_bc( dom.nodes[:(x==0 || x==1.0)], ux=0, uy=0)
set_bc( dom.nodes[:(y==0 || y==1.0)], ux=0, uy=0)

solve!(dom, nincs=20, verbose=true)

# boundary conditions
set_bc( dom.nodes[:(z==0.0)] , ux=0, uy=0, uz=0)
#set_bc( dom.nodes[:(z==0.5)] , uz=0.03)
set_bc( dom.faces[:(z==0.5)] , tz=2.)
set_bc( dom.nodes[:(x==0 || x==1.0)], ux=0, uy=0)
set_bc( dom.nodes[:(y==0 || y==1.0)], ux=0, uy=0)

solve!(dom, nincs=20, verbose=true)

# boundary conditions
set_bc( dom.nodes[:(z==0.0)] , ux=0, uy=0, uz=0)
set_bc( dom.nodes[:(z==0.5)] , uz=-0.035)
#set_bc( dom.nodes[:(x==0 || x==1.0)], ux=0, uy=0)
#set_bc( dom.nodes[:(y==0 || y==1.0)], ux=0, uy=0)

#solve!(dom, nincs=40, verbose=true)

#println(ip_dat[:ezz])
#println(ip_dat[:szz])
#save(dom, "dp.vtk")
save(ip_dat, "ip.dat")

#using Winston
#plot(ip_dat[:ezz], ip_dat[:szz])
