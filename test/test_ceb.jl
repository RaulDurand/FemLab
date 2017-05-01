using FemLab
using Base.Test

bl  = Block3D( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=1)
bli = BlockInset( [0.5 3 0.5; 0.5 6.0 0.5], curvetype="polyline")

mesh = Mesh(bl, bli, verbose=false)
dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=1.e4, nu=0.) )
set_state(dom.elems[:solids], sig=[-100, -100, -100, 0., 0., 0.])
#set_mat(dom.elems[:lines ], Truss(E=1.e8, A=0.005) )
set_mat(dom.elems[:lines ], Truss(E=1.e10, A=0.005) )

phi = 30*pi/180; dm=0.15138; c=20.0
#set_mat(dom.elems[:joints], CEBJoint1D(tau=[0. 0; 0.001 50; 0.0015 50; 0.003 5; 1000 5], kn=1.e5, dm=dm) )
set_mat(dom.elems[:joints1D], CEBJoint1D(ks=50000, TauR=10, s1=0.001, s2=0.0011, s3=0.002, A=0.005))

bar_nodes = get_nodes( dom.elems[:lines] )
bar_nodes = sort(bar_nodes, :y)
hook_node = bar_nodes[end]
solid_nodes = get_nodes(dom.elems[:solids])

tnode     = bar_nodes[end]
tab_tnode = NodeTracker(tnode)

tjoint     = dom.elems[:joints1D][end]
tab_tjoint = IpTracker(tjoint)

set_trackers(dom, tab_tnode, tab_tjoint)

## CEB test
scheme = "ME" 
disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=-0.0015)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=10, scheme=scheme, verbose=true)

disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=+0.0005)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=5, scheme=scheme, verbose=true)

disp_bc  = NodeBC( solid_nodes, ux=0, uy=0, uz=0)
force_bc = NodeBC( hook_node, uy=-0.0045)
set_bc(dom, disp_bc, force_bc)
@test solve!(dom, nincs=10, scheme=scheme, verbose=true)

#save(tab_tnode , "tab_tnode.dat")
#save(tab_tjoint, "tab_tjoint.dat")
