using FemLab
using Base.Test

bl  = Block3D( [0 0 0; 1. 1. 1.], nx=1, ny=1, nz=1, shape=HEX8)

# mesh generation
msh = Mesh(bl, verbose=false)

dom = Domain(msh) 
set_mat(dom.elems[:solids], Mazars(E=30000, nu=0.2, eps0=1.e-4, At=0.9, Bt=5000., Ac=1.0, Bc=1500.0))

# Tracking
bc1 = NodeBC( :(z==0), ux=0, uy=0, uz=0 )
bc2 = FaceBC( :(z==1), uz=-10e-3)
set_bc(dom, bc1, bc2)

@test solve!(dom, auto=true, nincs=20, scheme="FE", maxits=2, precision=0.01, verbose=true, saveincs=false)


