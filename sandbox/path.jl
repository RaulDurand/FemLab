#using FemLab

abstract IpData
abstract AbsJoint

include("mcjoint.jl")

mat = MCJoint(E=27e6, nu=0.2, ft=2.4e3, mu=1.4, alfa=0.4, beta=0.0, wc=1.7e-4, ws=1.85e-5 )
ipd = MCJointIpData()
ipd.h = 0.02
ipd.σ = [2400., 0, 0]

Δwt = [0.00001, 0.00002, 0]
nincs = 5

Δw = Δwt/nincs
for i=1:nincs
    print_with_color(:cyan, "Incremento $i:\n")
    Δσ = stress_update(mat, ipd, Δw)
end



