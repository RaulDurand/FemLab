using FemLab
using Base.Test


# 2D Truss

coord = [ 0 0; 1 0 ]
conn  = [ 1 2 ]

blt = BlockTruss(coord, conn)
mesh = Mesh(blt, verbose=false)

dom = Domain(mesh)

set_mat(dom.elems, Beam(E=10, A=1, I=1) )

bc1 = NodeBC( :(x==0 && y==0), ux=0, uy=0, rz=0)
bc2 = NodeBC( :(x==1 && y==0), fy=-10.)

set_bc(dom, bc1, bc2)

@test solve!(dom, verbose=true)

println(dom.nodes)
for node in dom.nodes
    println(node.dofs)
end
