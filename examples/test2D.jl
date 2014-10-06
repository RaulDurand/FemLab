using FemLab

bl = Block2D( [0 0; 1 1], nx=100, ny=100)
mesh = generate_mesh(bl)
#save(mesh, "mesh.vtk")

println(100)
dom = Domain(mesh)
println(200)
set_mat(dom.elems[:solids], ElasticSolid(E=1.0, nu=0.2) )
println(300)

set_bc( dom.nodes[:(x==0)] , ux=0, uy=0)
set_bc( dom.nodes[:(y==1)] , fy=-10.)



println(400)

solve!(dom, nincs=1, verbose=true)

