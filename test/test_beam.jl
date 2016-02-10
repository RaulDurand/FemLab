using FemLab
using FactCheck

verbose = true

bl  = Block3D( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=1)
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2] , curvetype="polyline")
bl2 = move( copy(bl1), x=0.6)

mesh = generate_mesh(bl, bl1, bl2, verbose=verbose)
save(mesh, "mesh.vtk")

dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=1.e4, nu=0.) )
set_mat(dom.elems[:lines ], Truss(E=1.e8, A=0.005) )
set_mat(dom.elems[:joints], Joint1D(ks=1.e5, kn=1.e5, A=0.005) )

bc1 = NodeBC( :(y==0 && z==0), uy=0, uz=0)
bc2 = NodeBC( :(y==6 && z==0), uz=0)
bc3 = FaceBC( :(z==1), tz=-10 )

set_bc(dom, bc1, bc2, bc3)

solve!(dom, nincs=1, verbose=verbose)
#save(dom, "viga.vtk")

facts("\nTest Reinforced Beam") do
    @fact 1 --> 1
end
