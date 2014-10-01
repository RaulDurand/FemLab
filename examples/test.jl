using FemLab

#bl1 = Block2D( [0 0; 1 1], nx=100, ny=100)
#bl2 = Block2D( [1 0; 2 1], nx=100, ny=100)
#bl3 = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10)

#blt = BlockTruss([0 0; 2 0; 1 1.], [1 2; 2 3; 1 3])

coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]

blt  = BlockTruss(coord, conn)
mesh = generate_mesh(blt)
save(mesh, "mesh.vtk")


#mesh = generate_mesh([bl1, bl2])
#mesh = generate_mesh(bl3)

dom = Domain(mesh)
#set_mat(dom.elems[:solids], MecElasticSolid(E=1.0, nu=0.2) )



set_mat(dom.elems, Truss(E=1.0, A=1.0) )

set_bc( dom.nodes[:(x==0 && y==0)] , ux=0, uy=0)
set_bc( dom.nodes[:(x==0 && y==9)] , ux=0, uy=0)
set_bc( dom.nodes[:(x==9 && y==0)] , fy=-450.)
set_bc( dom.nodes[:(x==18&& y==0)] , fy=-450.)


#set_bc(dom.nodes[:(x==2 && y==0)], uy=0)
#set_bc(dom.nodes[:(y==1)], fy=-10)
#set_bc(dom.nodes[:(z==0)], ux=0, uy=0, uz=0)
#set_bc(dom.nodes[:(z==1)], fz=-10)

#for n in dom.nodes[:(x==2)]
    #println("Node:")
    #for d in n.dofs
        #println(d, "\n")
    #end
#end
#@time begin 
#for i=1:100
    #println("ger. ", i)
    #for j=1:91
        solve!(dom, nincs=1, verbose=true)
    #end
#end
#end
#save(domain, "output1.vtk")

#A = [1,2,3]
#B = [1,2,3]

#C = 2*(A+B)

#for i in 1:n
    #C[i] = 2*A[i] + 2* B[i]
#end

readline()
