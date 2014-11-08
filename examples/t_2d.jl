using FemLab

# Mesh generation
bl = Block2D( [0 0; 1 1], nx=10, ny=30)
mesh = generate_mesh(bl)

save(mesh, "mesh.vtk")

# Finite element analysis
dom = Domain(mesh)

set_mat(dom.elems, ElasticSolid(E=1.0, nu=0.2) )

set_bc( dom.nodes[:(x==0)] , ux=0, uy=0)
set_bc( dom.nodes[:(y==1)] , fy=-10.)
set_bc( dom.faces[:(y==1)] , ty=-10.)

solve!(dom, nincs=1, verbose=true)
save(dom, "output.vtk")

