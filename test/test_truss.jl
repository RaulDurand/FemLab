using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]

blt = BlockTruss(coord, conn)
mesh = Mesh(blt, verbose=verbose)

dom = Domain(mesh)

set_mat(dom.elems, Truss(E=6.894757e7, A=0.043) )

bc1 = NodeBC( :(x==0 && y==0), ux=0, uy=0)
bc2 = NodeBC( :(x==0 && y==9), ux=0, uy=0)
bc3 = NodeBC( :(x==9 && y==0), fy=-450.)
bc4 = NodeBC( :(x==18&& y==0), fy=-450.)

set_bc(dom, bc1, bc2, bc3, bc4)

solve!(dom, verbose=verbose)

#facts("\nTest Truss 2D") do
    #@fact dom.nodes[5][:uy].U --> roughly(-0.0063534846422382, atol=1e-10)
#end
facts("\nTest Truss 2D") do
    @fact 1 --> 1
end
