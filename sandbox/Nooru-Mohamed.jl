using FemLab

# Material properties
E_c = 30000e3
f_t = 3.30e3

# Displacements
ws = 0.05*0.001
wt = 0.16*0.001

# Mesh
mesh = Mesh("Nooru-Mohamed.vtk")
#mesh = move(mesh, dx=-0.01)
split!(mesh)

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=E_c, nu=0.2) )
set_mat(dom.elems[:joints], MCJoint( E=E_c, nu=0.2, ft=f_t, mu=1.4, alfa=1.0,wc=1.70e-4*1, ws=1.85e-5*1))

#set_mat( [dom.elems[:solids][:(x<=0)]; dom.elems[:solids][:(y>=0.2)]], ElasticSolid(E=E_c*100, nu=0.0) )
#set_mat( [dom.elems[:joints][:(x<=0)]; dom.elems[:joints][:(y>=0.2)]], Joint(ks=1e9, kn=1e9) )

face_dat = FacesTracker( dom.faces[:(y==0.2)] )
set_trackers(dom, face_dat)


# Boundary conditions
bc1 = FaceBC(:(y==0), ux=0, uy=0, uz=0)
bc2 = FaceBC(:(x==0.2 && y<0.1), ux=0, uy=0, uz=0)

bc3 = FaceBC(:(y==0.2), uy=0, uz=0)
bc4 = FaceBC(:(x==0.0 && y>0.1), tx=10.0)

set_bc(dom, bc1, bc2, bc3, bc4)

# Solution
solve!(dom, auto=true, scheme="ME", precision=0.05, autosave=true)

# Boundary conditions
bc1 = FaceBC(:(y==0), ux=0, uy=0, uz=0)
bc2 = FaceBC(:(x==0.2 && y<0.1), ux=0, uy=0, uz=0)

bc3 = FaceBC(:(y==0.2), ux=0, uy=wt, uz=0)
bc4 = FaceBC(:(x==0.0 && y>0.1), ux=0, uy=wt)

set_bc(dom, bc1, bc2, bc3, bc4)

# Solution
solve!(dom, auto=true, nincs=20, scheme="ME", precision=0.05, savesteps=true)

save(face_dat, "face.dat")
