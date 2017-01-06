using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

bl  = Block3D( [0 0 0; 0.2 0.1 0.1], nx=2, ny=1, nz=1, shape=HEX8)

# mesh generation
msh = Mesh(bl, verbose=true)
split!(msh)

dom = Domain(msh)

E = 27.e6

set_mat(dom.elems[:solids], ElasticSolid(E=E, nu=0.2))
set_mat(dom.elems[:joints], MCJoint(E=E, nu=0.2, ft=2.4e3, mu=1.4, alfa=1.0, beta=1.0, wc=1.7e-4, ws=1.85e-5  ) )

# Tracking
midjoint = dom.elems[:joints][:(x==0.1)][1]
midjointdat = IpTracker(midjoint)
set_trackers(dom, midjointdat)

# Boundary conditions
bc1 = FaceBC( :(x==0), ux=0, uy=0, uz=0 )
#bc2 = FaceBC( :(x>0 ), uy=0, uz=0 )
bc3 = FaceBC( :(x==0.2), ux=2.0*1.7e-4)

set_bc(dom, bc1, bc3)
autosave = false

solve!(dom, auto=true, nincs=20, scheme="ME", precision=1e-3, autosave=autosave, verbose=verbose)

verbose && save(dom, "joint.vtk")


using PyPlot

plot(midjointdat.table[:w1], midjointdat.table[:s1], marker="^")
show()
#plot(midjointdat.table[:w1], midjointdat.table[:s2], marker="o")
#show()
#plot(midjointdat.table[:w1], midjointdat.table[:s3], marker="o")
#show()


facts("\nTest Joints") do
    @fact 1 --> 1
end

