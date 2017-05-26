using FemLab
using Base.Test

# mesh 
bl = Block3D( [0 0 0; 1 1 0.5], nx=2, ny=2, nz=2)
mesh = Mesh(bl, verbose=false)

# fem domain
dom = Domain(mesh)
set_mat(dom.elems, DruckerPrager(E=100., nu=0.25, alpha=0.05, kappa=0.1) )

ip_dat = IpTracker(dom.elems[1])
set_trackers(dom, ip_dat)

# boundary conditions
base  = NodeBC( :(z==0.0) , ux=0, uy=0, uz=0)
top   = NodeBC( :(z==0.5) , uz=-0.033)
sidex = NodeBC( :(x==0 || x==1.0), ux=0, uy=0)
sidey = NodeBC( :(y==0 || y==1.0), ux=0, uy=0)

set_bc(dom, base, top, sidex, sidey)
@test solve!(dom, auto=true, nincs=10, precision=1e-2, verbose=true)

# boundary conditions
top   = NodeBC( :(z==0.5) , uz=+0.008)
set_bc(dom, base, top, sidex, sidey)

@test solve!(dom, auto=true, nincs=10, verbose=true)

