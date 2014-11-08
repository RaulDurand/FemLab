using FemLab

function getline(file)
    line = ""
    while length(line)==0
        line = strip(readline(file))
    end
    return line
end

function parseline(line)
    vals = Any[]
    svals = split(line)
    for s in svals
        push!(vals, parse(s))
    end
    return vals
end

mesh = Mesh()
file = open("t_william.dat","r") 
getline(file)
getline(file)
getline(file)
nnodes, nelems, nsolids = parseline(getline(file))
@show nnodes
@show nelems
getline(file)
getline(file)
for i=1:nnodes
    n,x,y,z = parseline(getline(file))
    push!(mesh.points, Point(x,y,z) )
    mesh.points[i].id = i
end

for i=1:nelems
    con = parseline(getline(file))
    con = int(con[3:end])
    if i<=nsolids
        cell = Cell(HEX8, mesh.points[con])
    else
        cell = Cell(LIN2, mesh.points[con])
    end
    cell.id = i
    push!(mesh.cells, cell)
end

close(f)

save(mesh, "wmesh.vtk")
