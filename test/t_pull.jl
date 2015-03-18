using FemLab

bl  = Block3D( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=1)
bli = BlockInset( [0.5 2 0.5; 0.5 6.0 0.5], curvetype="polyline")

mesh = generate_mesh(bl, bli, verbose=false)
#save(mesh, "mesh.vtk")

dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=1.e4, nu=0.) )
set_state(dom.elems[:solids], sig=[-100, -100, -100, 0., 0., 0.])
set_mat(dom.elems[:lines ], Truss(E=1.e7, A=0.005) )
phi = 30*pi/180; dm=0.15138; c=20.0
set_mat(dom.elems[:joints], MCJoint1D(ks=1.e5, kn=1.e5, dm=dm, c=c, phi=phi) )

bar_nodes = get_nodes(dom.elems[:lines])
bar_nodes = sort(bar_nodes, :y)
hook_node = bar_nodes[end]

# Defining load levels
load_levels = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 0.9999]
load_incs   = diff(load_levels)
nstages = length(load_incs)
tload   = (pi*dm*4.)*(c+100*tan(phi))

solid_nodes = get_nodes(dom.elems[:solids])

# loop along stages
for i=1:nstages
    set_bc(solid_nodes, ux=0, uy=0, uz=0)
    set_bc(hook_node, fy = tload*load_incs[i])

    solve!(dom, nincs=1)
    #save(dom, "output$i.vtk")
end
