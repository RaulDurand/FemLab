using FemLab
using FactCheck

bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD8)
#bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD4)
mesh = generate_mesh(bl, verbose=false)

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=1.0, nu=0.2) )

set_bc( dom.nodes[:(y==0)] , ux=0, uy=0)
set_bc( dom.faces[:(y==1)] , ty=-10.)

solve!(dom, nincs=1, verbose=false)
#save(dom, "out8.vtk")

facts("\nTest 2D:") do
    @fact dom.nodes[end].dofdict[:uy].U --> roughly(-9.55, atol=6e-3)
end

