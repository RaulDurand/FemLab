using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

bl  = Block3D( [0 0 0; 3.0 0.3 0.4], nx=20, ny=1, nz=6, shape=HEX8)
bl1 = BlockInset( [0.05 0.05 0.05; 2.95 0.05 0.05], curvetype="polyline")
bl2 = move(copy(bl1), y=0.2)

# mesh generation
mesh = Mesh(bl, bl1, bl2, verbose=true)
split!(mesh)

dom = Domain(mesh)

#set_mat(dom.elems[:solids], ElasticSolid(E=27.e6, nu=0.2))
#set_mat(dom.elems[:joints], MCJoint(E=27e6, nu=0.2, ft=2.4e3, mu=1.4, alfa=0.4, beta=0.0, wc=1.7e-4, ws=1.85e-5  ) )
set_mat(dom.elems[:solids], ElasticSolid(E=38000e3, nu=0.2))
set_mat(dom.elems[:joints], MCJoint( E=38000e3, nu=0.2, ft=3.7e3, mu=1.4, alfa=1.0, beta=0.5 ,wc=1.17e-4, ws=2.28e-5))  

set_mat(dom.elems[:lines ], Truss(E=2.1e8, A=0.002) )
set_mat(dom.elems[:joints1D], Joint1D(ks=1.e8, kn=1.e8, A=0.002) )

# Tracking
#midjoint = dom.elems[:joints][:(x==0.1)][1]
#midjointdat = IpTracker(midjoint)
#set_trackers(dom, midjointdat)

# Boundary conditions
bc1 = NodeBC( :(x==0 && y==0 && z==0), uy=0 )
bc2 = NodeBC( :(x==0 && z==0), ux=0, uz=0 )
bc3 = NodeBC( :(x==3 && z==0), uz=0 )
bc4 = NodeBC( :(x==1.5 && z==0.4), uz=-0.01)

set_bc(dom, bc1, bc2, bc3, bc4)
#solve!(dom, nincs=20, precision=1e-1, autosave=verbose, verbose=verbose)
#solve_legacy!(dom, nincs=20, precision=1e-1, autosave=verbose, verbose=verbose)
solveFE!(dom, nincs=500, autosave=verbose, verbose=verbose)

verbose && save(dom, "joint.vtk")

facts("\nTest Joints") do
    @fact 1 --> 1
end
