using FemLab

bl = Block2D( [0 0; 1 1], nx=10, ny=10)
mesh = generate_mesh(bl)
#save(mesh, "mesh.vtk")

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=1.0, nu=0.2) )

set_bc( dom.nodes[:(x==0)] , ux=0, uy=0)
set_bc( dom.nodes[:(y==1)] , fy=-10.)


solve!(dom, nincs=1, verbose=true)

