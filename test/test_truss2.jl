using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

coord = [ 0. 0.; 2. 0.; 1. 1. ]
conn  = [ 1 2; 2 3; 1 3]

blt = BlockTruss(coord, conn)
mesh = Mesh(blt, verbose=verbose)

dom = Domain(mesh)

set_mat(dom.elems, Truss(E=1., A=1.) )

bc1 = NodeBC( :(x==0 && y==0), ux=0, uy=0)
bc2 = NodeBC( :(x==2 && y==0), uy=0)
bc3 = NodeBC( :(x==1. && y==1.), fy=-10.)

set_bc(dom, bc1, bc2, bc3)

tnode  = NodeTracker(dom.nodes[3])
tnodes = NodesTracker(dom.nodes)
set_trackers(dom, tnode, tnodes)

solve!(dom, nincs=10, verbose=verbose)

verbose && save(tnode, "tab.dat")
verbose && save(tnodes, "tnodes.dat")

facts("\nTest Truss 2D") do
    @fact 1 --> 1
end
