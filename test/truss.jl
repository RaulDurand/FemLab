using FemLab
using Base.Test


# 2D Truss

coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]

blt = BlockTruss(coord, conn)
mesh = Mesh(blt, verbose=false)

dom = Domain(mesh)

set_mat(dom.elems, Truss(E=6.894757e7, A=0.043) )

bc1 = NodeBC( :(x==0 && y==0), ux=0, uy=0)
bc2 = NodeBC( :(x==0 && y==9), ux=0, uy=0)
bc3 = NodeBC( :(x==9 && y==0), fy=-450.)
bc4 = NodeBC( :(x==18&& y==0), fy=-450.)

set_bc(dom, bc1, bc2, bc3, bc4)

@test solve!(dom, verbose=true)


# 3D Truss

coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]   # matriz de coordenadas
conn  = [ 1 3; 1 2;2 3]  # matriz de conectividades

blt = BlockTruss(coord, conn)
mesh = Mesh(blt, verbose=false)

dom = Domain(mesh)

set_mat(dom.elems, Truss(E=2.0e8, A=0.002) )

bc1 = NodeBC(:(x==0 && y==0 && z==0), ux=0)
bc2 = NodeBC(:(x==0 && y==1 && z==0), ux=0, uy=0, uz=0)
bc3 = NodeBC(:(x==0 && y==1 && z==1), ux=0, uy=0)
bc4 = NodeBC(:(x==0 && y==0), fz=-50.)

set_bc(dom, bc1, bc2, bc3, bc4)

@test solve!(dom, verbose=false)

