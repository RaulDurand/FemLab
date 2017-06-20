using FemLab
using Base.Test

bl = Block3D( [0 0 0; 1 1 1], nx=4, ny=4, nz=4, shape=HEX8)
mesh = Mesh(bl, verbose=false)

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=100.0, nu=0.2) )

node_mon1 = NodeMonitor( :(x==1 && y==1 && z==1) )
node_mon2 = NodeMonitor( :(x==1 && y==0 ) )
node_mon3 = NodeMonitor( dom.nodes[1] )

nodes_mon1 = NodesMonitor( dom.nodes )
nodes_mon2 = NodesMonitor( :(x==1 && y==1 ) )

faces_mon1 = FacesMonitor( dom.faces[:(z==1)] )
faces_mon2 = FacesMonitor( :(z==1) )
edges_mon1 = EdgesMonitor( dom.edges[:(z==1)] )
edges_mon2 = EdgesMonitor( :(z==1) )

ip_mon1   = IpMonitor( dom.elems[1].ips[1] )
#ip_mon2   = IpMonitor( dom.elems[1] )
ip_mon3   = IpMonitor( :(x>0.5 && y>0.5 && z>0.5) )

ips_mon1   = IpsMonitor( dom.elems[1].ips )
#ips_mon2   = IpsMonitor( dom.elems )
ips_mon3   = IpsMonitor( :(x>0.5 && y>0.5) )

set_monitors(dom, node_mon1, node_mon2, node_mon3, nodes_mon1, nodes_mon2, faces_mon1, faces_mon2,
             edges_mon1, edges_mon2, ip_mon1, ip_mon3, ips_mon1, ips_mon3)

base_bc = NodeBC( :(z==0), ux=0, uy=0, uz=0 )
lat_bc  = NodeBC( :(x==0 || x==1), ux=0)
top_bc  = FaceBC( :(z==1), tz=-10. )
set_bc(dom, base_bc, top_bc, lat_bc)

@test solve!(dom, nincs=1, verbose=true)
println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].U)
