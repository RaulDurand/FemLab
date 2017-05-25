using FemLab
using Base.Test

bl  = Block3D( [0 0 0; 1.0 6.0 1.0], nx=1, ny=60, nz=1)
bli = BlockInset( [0.5 3 0.5; 0.5 6.0 0.5], curvetype="polyline")

msh = Mesh(bl, bli, verbose=false)
dom = Domain(msh)

set_mat(dom.elems[:solids], ElasticSolid(E=24e3, nu=0.2) )
set_state(dom.elems[:solids], sig=[-100, -100, -100, 0., 0., 0.])
#set_mat(dom.elems[:lines ], Truss(E=1.e8, A=0.005) )
set_mat(dom.elems[:lines ], Truss(E=200e5, A=0.00011) )

phi = 30*pi/180; dm=0.15138; c=20.0
#set_mat(dom.elems[:joints], CEBJoint1D(tau=[0. 0; 0.001 50; 0.0015 50; 0.003 5; 1000 5], kn=1.e5, dm=dm) )
set_mat(dom.elems[:joints1D], CEBJoint1D(TauM=12, TauR=3, s1=0.001, s2=0.0011, s3=0.006, alpha=0.5, beta=0.5, ks=(12/0.001)*5,  kn=50000, A=0.005))

bar_nodes = dom.elems[:lines][:nodes]
bar_nodes = sort(bar_nodes, :y)
hook_node = bar_nodes[end]
solid_nodes = dom.elems[:solids][:nodes]


joint_ips = dom.elems[:joints1D][:ips]
joint_ips = sort(joint_ips, :y)
mon_joint_ips = IpsMonitor(joint_ips)


tnode     = bar_nodes[end]
tab_tnode = NodeTracker(tnode)

tjoint     = dom.elems[:joints1D][end]
mon_tjoint = IpMonitor(tjoint)

set_monitors(dom, tab_tnode, mon_tjoint, mon_joint_ips)

## CEB test
scheme = "FE" 
disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=+0.0003)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=40, auto=true, scheme=scheme, verbose=true)

disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=-0.0001)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=20, scheme=scheme, verbose=true)

disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=0.0006)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=40, scheme=scheme, verbose=true)

disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=-0.0005)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=20, scheme=scheme, verbose=true)

disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=0.007)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=60, scheme=scheme, verbose=true)

#using PyPlot
#tab = mon_tjoint.table
#save(mon_tjoint, "tab.dat")
#plot(tab[:ur], tab[:tau], marker="o", color="blue")
#show()
#book = mon_joint_ips.book
#save(mon_joint_ips, "book.dat")
#for tab in book.tables
    #plot(tab[:y], tab[:tau], marker="o")
#end
#show()
