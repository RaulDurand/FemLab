using FemLab
using Base.Test

# Mesh generation
bl  = Block3D( [0 0 0; 1.0 6.0 1.0], nx=1, ny=10, nz=1)
bl1 = BlockInset( [0.2 0.2 0.2; 0.2 5.8 0.2] , curvetype="polyline")
bl2 = move( copy(bl1), x=0.6)

mesh = Mesh(bl, bl1, bl2, verbose=true)

generate_embedded_cells!(mesh)  ## Deletes joint1D cells and sets up line cells as embedded cells

# FEM analysis
dom = Domain(mesh)

set_mat(dom.elems[:solids], ElasticSolid(E=1.e4, nu=0.) )
set_mat(dom.elems[:embedded], EmbPPTruss(E=1.e8, A=0.005, sig_y=500e3) )

bc1 = NodeBC( :(y==0 && z==0), ux=0, uy=0, uz=0)
bc2 = NodeBC( :(y==6 && z==0), uz=0)
bc3 = FaceBC( :(z==1), tz=-1000 )

set_bc(dom, bc1, bc2, bc3)

@test solve!(dom, nincs=20, verbose=true, saveincs=false)

