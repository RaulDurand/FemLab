using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

coord = [ 0. 0.; 1. 0. ]
conn  = [ 1 2 ]

blt  = BlockTruss(coord, conn)
mesh = Mesh(blt, verbose=verbose)

dom = Domain(mesh)

set_mat(dom.elems, PPTruss(E=210e6, A=0.01, sig_y=500e3) )
#set_mat(dom.elems, PPTruss(E=210e6, A=0.01, sig_y=500e3, H=1000) )

bc1 = NodeBC( :(x==0 && y==0), ux=0, uy=0)
bc2 = NodeBC( :(x==1 && y==0), uy=0)
bc3 = NodeBC( :(x==1 && y==0), ux=0.003)
#bc3 = NodeBC( :(x==1 && y==0), fx=5010.000)

set_bc(dom, bc1, bc2, bc3)

#tnode  = NodeTracker(dom.nodes[3])
#tnodes = NodesTracker(dom.nodes)
#set_trackers(dom, tnode, tnodes)

solve!(dom, nincs=10, verbose=verbose)

#verbose && save(tnode, "tab.dat")
#verbose && save(tnodes, "tnodes.dat")

facts("\nTest Truss 2D") do
    @fact 1 --> 1
end
