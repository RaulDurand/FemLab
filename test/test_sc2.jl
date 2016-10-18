using FemLab
using FactCheck

verbose = isdefined(:verbose) ? verbose : true

bl  = Block3D( [0 0 0; 0.2 0.8 0.2], nx=1, ny=1, nz=1)

mesh = Mesh(bl, verbose=verbose)
dom = Domain(mesh)

set_mat(dom.elems[:solids], SmearedCrack(E=20e6, nu=0.20, fc=48e3, ft=2.4e3, Gf=0.113, xi1=0.4, xi2=0.8, al1=0.6, al2=0.2, p1=2, h=0.01) )

node_dat = NodeTracker( dom.nodes[:(x==0.2 && y==0.8 && z==0.2)] )
set_trackers(dom, node_dat)


bc1 = NodeBC( :(y==0), ux=0, uy=0, uz=0)
bc2 = NodeBC( :(y==0.8), uy=0.01)

set_bc(dom, bc1, bc2)

solve!(dom, nincs=20, verbose=verbose, autosave=true)

using PyPlot
plot(node_dat.table[:uy], node_dat.table[:fy], marker="o")
show()
