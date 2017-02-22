using FemLab

bl  = Block3D( [0 0 0; 1. 1. 1.], nx=1, ny=1, nz=5, shape=HEX8)

# mesh generation
msh = Mesh(bl, verbose=true)


dom = Domain(msh) # Units in MPa
set_mat(dom.elems[:solids], Mazars(E=30000, nu=0.2, eps0=1.e-4, At=0.9, Bt=5000., Ac=1.0, Bc=1500.0))

# Tracking
ip1 = dom.elems[end][:ips][1]
ip1dat = IpTracker(ip1, "ip1.dat")
set_trackers(dom, ip1dat)


bc1 = NodeBC( :(z==0), ux=0, uy=0, uz=0 )
bc2 = FaceBC( :(z==1), uz=-10e-3)
set_bc(dom, bc1, bc2)

try solve!(dom, auto=true, nincs=40, scheme="FE", maxits=2, precision=0.01, autosave=true)
catch err @show err end


#bc1 = NodeBC( :(z==0), ux=0, uy=0, uz=0 )
#bc2 = FaceBC( :(z==1), uz=+3e-3)
#set_bc(dom, bc1, bc2)

#try 
    #solve!(dom, auto=true, nincs=100, scheme="FE", maxits=2, precision=0.01, autosave=true)
#catch err @show err end


#exit()

using PyPlot

tab = ip1dat.table
plot(tab[:ezz], tab[:szz], marker="^")
plot(tab[:ezz], tab[:dam], marker="o")
show()
#plot(tab[:eq], tab[:dam], marker="o")
#show()
