using FemLab
using Base.Test

bl = Block2D( [0 0; 1 1], nx=8, ny=8, shape=QUAD8)
mesh = Mesh(bl, verbose=false)

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=1.0, nu=0.2) )

base_bc = NodeBC( :(y==0), ux=0, uy=0 )
top_bc  = FaceBC( :(y==1), ty=-10. )

set_bc(dom, base_bc, top_bc)

@test solve!(dom, nincs=1, verbose=true)

