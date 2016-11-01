using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

bl  = Block3D( [0 0 0; 0.2 0.1 0.1], nx=2, ny=1, nz=1, shape=HEX8)

# mesh generation
mesh = Mesh(bl, verbose=true)
split!(mesh)

dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=27.e6, nu=0.2))
set_mat(dom.elems[:joints], MCJoint(E=27e6, nu=0.2, ft=2.4e3, mu=1.4, alfa=0.4, beta=1.0, wc=1.7e-4, ws=1.85e-5  ) )

# Tracking
midjoint = dom.elems[:joints][:(x==0.1)][1]
midjointdat = IpTracker(midjoint)
set_trackers(dom, midjointdat)

# Boundary conditions
bc1 = FaceBC( :(x==0), ux=0, uy=0, uz=0 )
bc2 = FaceBC( :(x==0.2), uy=0.9*1.7e-4)

set_bc(dom, bc1, bc2)
solve!(dom, nincs=20, precision=1e-4, autosave=verbose, verbose=verbose)

verbose && save(dom, "joint.vtk")

facts("\nTest Joints") do
    @fact 1 --> 1
end
