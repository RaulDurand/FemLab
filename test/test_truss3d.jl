using FemLab
using FactCheck

verbose=true

coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]   # matriz de coordenadas
conn  = [ 1 3; 1 2;2 3]  # matriz de conectividades

blt = BlockTruss(coord, conn)
mesh = generate_mesh(blt, verbose=verbose)

dom = Domain(mesh)

set_mat(dom.elems, Truss(E=2.0e8, A=0.002) )

bc1 = NodeBC(:(x==0 && y==0 && z==0), ux=0)
bc2 = NodeBC(:(x==0 && y==1 && z==0), ux=0, uy=0, uz=0)
bc3 = NodeBC(:(x==0 && y==1 && z==1), ux=0, uy=0)
bc4 = NodeBC(:(x==0 && y==0), fz=-50.)

set_bc(dom, bc1, bc2, bc3, bc4)

solve!(dom, verbose=verbose)

facts("\nTest Truss 3D") do
    @fact 1 --> 1
end
