using FemLab
using FactCheck

bl = Block2D( [0 0; 1 1], nx=10, ny=10)
mesh = generate_mesh(bl, verbose=false)

dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=1.0, nu=0.2) )

set_bc( dom.nodes[:(y==0)] , ux=0, uy=0)
set_bc( dom.nodes[:(y==1)] , fy=-10.)

solve!(dom, nincs=1, verbose=false)

facts("\nTest 2D:") do
    @fact dom.nodes[121].dofdict[:uy].U => roughly(-127.37914254439363, atol=1e-10)
end

