using FemLab
using FactCheck

#bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=FemMesh.QUAD8)
bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD8)
#bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD4)
mesh = generate_mesh(bl, verbose=false)

#mesh = load_mesh("file.vtk", format="vtk")

dom = Domain(mesh)
#dom = load_domain("input.json")
#dom = Domain("input.json")

set_mat(dom.elems[:solids], ElasticSolid(E=1.0, nu=0.2) )

set_bc( @get_nodes(dom, y==0), ux=0, uy=0)
set_bc( @get_faces(dom, y==1), ty=-10.)

#set_bc( dom.nodes[:(y==0)] , ux=0, uy=0)
#set_bc( dom.faces[:(y==1)] , ty=-10.)

solve!(dom, nincs=1, verbose=false)

facts("\nTest 2D:") do
    @fact dom.nodes[end].dofdict[:uy].U --> roughly(-9.55, atol=6e-3)
end

