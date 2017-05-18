using FemLab
using Base.Test

for shape in (TRI3, TRI6, QUAD4, QUAD8, QUAD9)
    bl = Block2D( [0 0; 1 1], nx=8, ny=8, shape=shape)
    mesh = Mesh(bl, verbose=false)

    dom = Domain(mesh)
    set_mat(dom.elems[:solids], ElasticSolid(E=100.0, nu=0.2) )

    base_bc = NodeBC( :(y==0), ux=0, uy=0 )
    top_bc  = FaceBC( :(y==1), ty=-10. )

    set_bc(dom, base_bc, top_bc)

    @test solve!(dom, nincs=1, verbose=true)
    println("  uy = ", dom.nodes[:(y==1)][1].dofdict[:uy].U)
end

for shape in (TET4, TET10, HEX8, HEX20)
    bl = Block3D( [0 0 0; 1 1 1], nx=4, ny=4, nz=4, shape=shape)
    mesh = Mesh(bl, verbose=false)

    dom = Domain(mesh)
    set_mat(dom.elems[:solids], ElasticSolid(E=100.0, nu=0.2) )

    base_bc = NodeBC( :(z==0), ux=0, uy=0, uz=0 )
    lat_bc  = NodeBC( :(x==0 || x==1), ux=0)
    top_bc  = FaceBC( :(z==1), tz=-10. )

    set_bc(dom, base_bc, top_bc, lat_bc)

    @test solve!(dom, nincs=1, verbose=true)
    println("  uz = ", dom.nodes[:(z==1)][1].dofdict[:uz].U)
end
