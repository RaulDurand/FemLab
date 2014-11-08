using FemLab

bl  = Block3D( [0 0 0; 3. 0.25 0.42], shape=HEX20, nx=12, ny=3, nz=6)
bl1 = BlockInset( [0.05 0.03 0.02; 2.95 0.03 0.02], curvetype="polyline")
bl2 = move(copy(bl1), y=0.08)
bl3 = move(copy(bl1), y=0.14)
bl4 = move(copy(bl1), y=0.19)

#mesh = generate_mesh(bl)
mesh = generate_mesh(bl, bl1, bl2, bl3, bl4)
save(mesh, "mymesh.vtk")
#exit()

dom = Domain(mesh)
A = 0.0001

set_mat(dom.elems[:solids], DruckerPrager(E=10.e6, nu=0.2, alpha=0.42, kappa=6000., H=0.) )
set_mat(dom.elems[:lines ], PPTruss(E=2.1e8, A=A, sig_y=434.e3) )
set_mat(dom.elems[:joints], Joint1D(ks=1.e9, kn=1.e9, A=A) )

t_node9 = track(dom, dom.nodes[9])
t_elem9 = track(dom, dom.elems[9])

#w=-650 # load
w=-400 # load
for i=1:2
    set_bc(dom.nodes[:(x==0 && z==0)], ux=0, uy=0, uz=0)
    set_bc(dom.nodes[:(x==3 && z==0)], uy=0, uz=0)
    set_bc(dom.faces[:(z==0.42)], tz=0.5*w)

    solve!(dom, nincs=4)
    save(dom, "output$i.vtk")
end

save(t_node9, "t_node9.dat")
save(t_elem9, "t_elem9.dat")
