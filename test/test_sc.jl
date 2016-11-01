using FemLab
using FactCheck

verbose = isdefined(:verbose) ? verbose : true

bl  = Block3D( [0 0 0; 0.3 6.0 0.4], nx=1, ny=40, nz=6)
bl1 = BlockInset( [0.05 0.05 0.05; 0.05 5.95 0.05] , curvetype="polyline")
bl2 = move( copy(bl1), x=0.1)
bl3 = move( copy(bl1), x=0.2)

#mesh = Mesh(bl, bl1, bl2, bl3, verbose=verbose)
mesh = Mesh(bl, verbose=verbose)

dom = Domain(mesh)

#set_mat(dom.elems[:solids], ElasticSolid(E=1.e4, nu=0.) )
set_mat(dom.elems[:solids], SmearedCrack(E=20e6, nu=0.20, fc=48e3, ft=2.4e3, Gf=0.113, xi1=0.4, xi2=0.8, al1=0.6, al2=0.2, p1=2, h=0.01) )
#set_mat(dom.elems[:lines ], Truss(E=2.1e8, A=0.01) )
#set_mat(dom.elems[:joints1D], Joint1D(ks=1.e5, kn=1.e5, A=0.005) )

mid_node = dom.nodes[:(y==3.0 && z==0.4)][1]
node_dat = NodeTracker( mid_node )
set_trackers(dom, node_dat)

bc1 = NodeBC( :(y==0 && z==0), ux=0, uy=0, uz=0)
bc2 = NodeBC( :(y==6 && z==0), ux=0, uz=0)
#bc3 = FaceBC( :(z==0.4), tz=-30)
bc3 = NodeBC( :(y==3 && z==0.4), uz=-0.003)
#bc3 = NodeBC( :(y==3.0 && z==0.4), uz=-0.01)

set_bc(dom, bc1, bc2, bc3)

solve!(dom, nincs=10, verbose=verbose, autosave=true)

#save(dom, "")
if verbose
    using PyPlot
    #plot(-node_dat.table[:uz], 0:80, marker="o")
    plot(-node_dat.table[:uz], -node_dat.table[:fz], marker="o")
    show()
end
