using FemLab

# Material properties
E_c = 30000e3
f_t = 3.30e3

# Displacements
ws = 0.05*0.001
wt = 0.16*0.001

# Mesh
mesh = Mesh("Nooru-Mohamed.vtk")
mesh = Mesh(mesh)
split!(mesh)
#save(mesh, "mesh.vtk")

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=E_c, nu=0.2) )

set_mat(dom.elems[:joints], MCJoint( E=E_c, nu=0.2, ft=f_t, mu=1.4, alfa=1.0,wc=1.70e-4*1, ws=1.85e-5*1))

# Boundary conditions
bc1 = FaceBC(:(y==0), ux=0, uy=0, uz=0)
bc2 = FaceBC(:(x==0.2 && y<0.1), ux=0, uy=0, uz=0)

bc3 = FaceBC(:(y==0.2), uy=0, uz=0)
bc4 = FaceBC(:(x==0 && y>0.1), ux=ws)

set_bc(dom, bc1, bc2, bc3, bc4)

# Solution
solve!(dom, auto=true, scheme="ME", precision=1e-2, autosave=true)

# Boundary conditions
bc1 = FaceBC(:(y==0), ux=0, uy=0, uz=0)
bc2 = FaceBC(:(x==0.2 && y<0.1), ux=0, uy=0, uz=0)

bc3 = FaceBC(:(y==0.2), ux=0, uy=wt, uz=0)
bc4 = FaceBC(:(x==0 && y>0.1), ux=0, uy=wt)

set_bc(dom, bc1, bc2, bc3, bc4)

# Solution
solve!(dom, auto=true, nincs=20, scheme="ME", precision=1e-2, savesteps=true)

