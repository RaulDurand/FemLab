using FemLab
n = Node(rand(3))
l = repeat([n], inner=100)
@show l
#println(l)
