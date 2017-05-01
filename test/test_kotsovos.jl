using FemLab
using Base.Test

bl  = Block3D( [0 0 0; 0.3 6.0 0.4], nx=1, ny=20, nz=4)
bl1 = BlockInset( [0.05 0.05 0.05; 0.05 5.95 0.05] , curvetype="polyline")
bl2 = move( copy(bl1), x=0.1)
bl3 = move( copy(bl1), x=0.2)

mesh = Mesh(bl, bl1, bl2, bl3, verbose=false)

dom = Domain(mesh)

set_mat(dom.elems[:solids], Kotsovos(E=20e6, nu=0.25, beta=0.5, fc=20e3, ft=2e3) )
set_mat(dom.elems[:lines ], Truss(E=2.1e8, A=0.01) )
set_mat(dom.elems[:joints1D], Joint1D(ks=1.e5, kn=1.e5, A=0.005) )

mid_node = dom.nodes[:(y==3.0 && z==0.4)][1]
node_dat = NodeTracker( mid_node )
set_trackers(dom, node_dat)

bc1 = NodeBC( :(y==0 && z==0), ux=0, uy=0, uz=0)
bc2 = NodeBC( :(y==6 && z==0), ux=0, uz=0)
bc3 = NodeBC( :(y==3.0 && z==0.4), uz=-0.005)

set_bc(dom, bc1, bc2, bc3)

@test solve!(dom, nincs=40, precision=0.1, verbose=false)

#solve_legacy!(dom, nincs=50, precision=0.001, verbose=false, autosave=true)
#solveFE!(dom, nincs=50, verbose=false, autosave=true)

#using PyPlot
#plot(-node_dat.table[:uz], -node_dat.table[:fz], marker="o")
#show()
