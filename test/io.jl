using FemLab

#=
for shape in (TRI3, TRI6, QUAD4, QUAD8, QUAD9)
    bl = Block2D( [0 0; 1 1], nx=2, ny=2, shape=shape)
    mesh = Mesh(bl, verbose=false)

    dom = Domain(mesh)
    save(dom, "dom.vtk")
end
rm("dom.vtk")
=#

bl = Block2D( [0 0; 1 1], nx=4, ny=4, shape=QUAD9)
mesh = Mesh(bl, verbose=false)
dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=100.0, nu=0.2) )
bc1 = NodeBC(:(y==0), ux=0, uy=0)
bc2 = FaceBC(:(y==1), tz=2)

println(dom.nodes)
for node in dom.nodes[1:5]
    println(node.dofs)
end
println(dom.elems)
println(dom.elems[:ips])

println(dom)
println(bc1)
println(bc2)
