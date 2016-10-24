using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

bl  = Block3D( [0 0 0; 0.2 0.1 0.1], nx=2, ny=1, nz=1, shape=HEX8)

mesh = Mesh(bl, verbose=true)
split!(mesh)

dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=27.e6, nu=0.2))

#set_mat(dom.elems[:joints], Joint(ks=1.e4, kn=2.e4) )
set_mat(dom.elems[:joints], MCJoint(E=27e6, nu=0.2, ft=2.4e3, mu=1.4, alfa=5.0, beta=1.0, wc=1.7e-4, ws=1.85e-5  ) )

fixl = NodeBC( :(x==0), ux=0, uy=0, uz=0 )
#fixr = NodeBC( :(x==0.2), uy=0, uz=0 )
#load = FaceBC( :(z==1.0), tz=-150.0)
load = NodeBC( :(x==0.2), ux=0.0002)
set_bc(dom, load, fixl)

solve!(dom, nincs=10, precision=1e-4, autosave=verbose, verbose=verbose)

verbose && save(dom, "joint.vtk")

facts("\nTest Joints") do
    @fact 1 --> 1
end
