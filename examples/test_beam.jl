using FemLab

bl  = Block3D( [0 0 0; 0.4 4.0 0.4], nx=10, ny=50, nz=20)
bli = BlockInset( [0.2 0 0.1; 0.2 4.0 0.1], curvetype="polyline")

mesh = generate_mesh(bl, bli)
save(mesh, "outmesh.vtk")

dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=10.0, nu=0.2) )
set_mat(dom.elems[:lines ], PPTruss(E=1.9, A=0.1, sig_y=8.0) )
set_mat(dom.elems[:joints], Joint1D(ks=1.0, kn=0.1, dm=0.1) )

set_bc(dom.nodes[:(y==0 && z==0)], ux=0, uy=0, uz=0)
set_bc(dom.nodes[:(y==4 && z==0)], ux=0, uz=0)
set_bc(dom.faces[:(z==0.4)], tz=-10)

solve!(dom, nincs=1)
save(dom, "output1.vtk")

