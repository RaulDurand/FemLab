using FemLab
using Base.Test

bl  = Block2D( [0 0; 0.2 0.1], nx=2, ny=1, shape=QUAD4)

# mesh generation
msh = Mesh(bl, verbose=false)
split!(msh)

dom = Domain(msh, stress_state=:plane_stress, thickness=1.0)

E = 27.e6

set_mat(dom.elems[:solids], ElasticSolid(E=E, nu=0.2))
set_mat(dom.elems[:joints], MCJoint(E=E, nu=0.2, ft=2.4e3, mu=1.4, alpha=1.0, wc=1.7e-4, ws=1.85e-5 ) )

# Tracking
midjoint = dom.elems[:joints][:(x==0.1)][1][:ips][2]

# Boundary conditions
bc1 = FaceBC( :(x==0), ux=0, uy=0 )
bc3 = FaceBC( :(x==0.2), ux=2.0*1.7e-4)

set_bc(dom, bc1, bc3)

@test solve!(dom, auto=true, nincs=20, scheme="ME", maxits=3, precision=0.01, verbose=true)
