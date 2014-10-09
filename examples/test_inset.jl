using FemLab

bl  = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10)
bli = BlockInset( [0 0 0; 0.3 0.3 0.7; 0.7 0.7 0.3; 1 1 1], curvetype="Bezier")

mesh = generate_mesh(bl, bli)
save(mesh, "outmesh.vtk")
exit()

dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=20.0, nu=0.2) )
set_mat(dom.elems[:lines ], Truss(E=1.0, A=0.1) )
#set_mat(dom.elems[:lines ], PPTruss(E=20.0, A=0.1, sig_y=0.01) )
set_mat(dom.elems[:joints], Joint1D(ks=1.0, kn=0.1, dm=0.1) )

set_bc( dom.faces[:(z==0)] , ux=0, uy=0, uz=0)
set_bc( dom.faces[:(z==1)] , tz=-10)

#datae = track(dom, dom.elems[1])
#datan = track(dom, dom.nodes[1])
#datac = track(dom, dom.nodes[1:4])
#dataip = track(dom, dom.elems[1].ips[1])

solve!(dom, nincs=1)
save(dom, "output1.vtk")

