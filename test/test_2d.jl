using FemLab
using FactCheck

verbose = true

bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD8)
mesh = generate_mesh(bl, verbose=verbose)

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=1.0, nu=0.2) )

base_bc = NodeBC( :(y==0), ux=0, uy=0 )
top_bc  = FaceBC( :(y==1), ty=-10. )

set_bc(dom, base_bc, top_bc)

solve!(dom, nincs=1, verbose=verbose)

save(dom, "out.vtk")

