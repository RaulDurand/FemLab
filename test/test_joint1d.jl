using FemLab
using FactCheck
verbose = isdefined(:verbose) ? verbose : true

#bl  = Block3D( [0 0 0; 1.0 1.0 1.0], nx=20)
bl  = Block3D( [ -0.05 -0.05 0; 0.05 0.05 0.2 ], nx=1, ny=1, nz=2, shape=HEX8 )
#bli = BlockInset( [0.5 0.5 0.5; 1.0 0.5 0.5])
bli = BlockInset( [0.01 0.010 0.1; 0.01 0.01 0.2] )

msh = Mesh(bl, bli, verbose=verbose)

dom = Domain(msh)

set_mat(dom.elems[:solids], ElasticSolid(E=1.e4, nu=0.) )
set_mat(dom.elems[:lines ], Truss(E=1.e5, A=0.001) )
#set_mat(dom.elems[:joints1D], Joint1D(ks=1.e6, kn=1.e6, A=0.001) )
set_mat(dom.elems[:joints1D], CEBJoint1D(ks=50000, TauR=10, s1=0.001, s2=0.0011, s3=0.002, A=0.001))

tjoint   = IpTracker(dom.elems[:joints1D][end])
set_trackers(dom, tjoint)

disp_bc  = NodeBC(dom.elems[:solids][:nodes][:(z==0.0)], ux=0, uy=0, uz=0)
force_bc = NodeBC(dom.elems[:lines][:nodes][:(z==0.2)], uz= 0.01 )
set_bc(dom, disp_bc, force_bc)

solve!(dom, nincs=40, autosave=true, verbose=verbose)

verbose && save(dom, "dom.vtk")
verbose && save(tjoint, "tjoint.dat")

facts("\nTest Joint1D") do
    @fact 1 --> 1
end

#!verbose && exit()

#using PyPlot
#plot(tjoint.table[:ur], tjoint.table[:tau], marker="o")
#show()
