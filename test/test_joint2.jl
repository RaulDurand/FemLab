using FemLab

#Geometria da viga

bl = Block3D( [0 0 0; 0.3 0.8 0.3], shape=HEX8, nx=1, ny=2, nz=1)

# geracao malha
mesh = Mesh(bl)
split!(mesh) 
save(mesh, "mesh.vtk")

dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=38000e3,nu=0.2) )
set_mat(dom.elems[:joints], MCJoint( E=38000e3, nu=0.2, ft=3.7e3, mu=1.4, alfa=1, beta=1 ,wc=1.17e-4, ws=1.28e-5))  


# Acompanhamento de nós

midjoint = dom.elems[:joints][:(x==0.4)]

#midjointdat = IpTracker(midjoint)
#set_trackers(dom, midjointdat)

# Condicoes de contorno e solucao
bc1 = FaceBC(:(y==0), ux=0, uy=0, uz=0)
bc2 = FaceBC(:(y==0.8), uy=2.0*1.17e-4)
#bc2 = FaceBC(:(y==0.8), ux=20*234e-6)

set_bc(dom, bc1, bc2)
solve!(dom, nincs=40, autosave=true)

save(dom, "domain.vtk")
#save(nomeio_dat, "no719.dat")
#save(midjointdat, "joint.dat")

#using PyPlot

#plot(midjointdat.table[:w], midjointdat.table[:s_n], marker="o")
#show()
#plot(midjointdat.table[:s1], midjointdat.table[:tau1], marker="o")
#show()
#plot(midjointdat.table[:s2], midjointdat.table[:tau2], marker="o")
#show()
