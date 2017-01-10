using FemLab

# Geometria da viga

# Dimenções elementos da viga
Eviga    =  0.0800		   # espessura viga
hviga    =  0.0630		   # altura viga
Lviga    =  0.3500		   # largura viga
Lparcial =  0.2500		   # largura entre apoios
lfiss    =  0.0010         # espessura da fissura
hapoio   =  0.0025		   # altura apoio

#Propriedades do concreto
E_c     =  32000e3		   # modulo de elasticidade do concreto (kPa)
f_t     =  4.15e3		   # resistencia a tração concreto (kPa)
f_c     =  58.3e3		   # resistencia a compressão concreto (kPa)

#Propriedades apoio
E_s     =  218000e3		   # modulo de elasticidade do aço (kPa)

# Deslocamentos

w = 0.0009                    # deslocamento aplicado (m)


# Geracao malha

mesh = Mesh("roesler.vtk")
#mesh  = Mesh(bviga, genedges=true)

split!(mesh)


# Definicao dominio: materiais
dom = Domain(mesh)
set_mat(dom.elems[:solids], ElasticSolid(E=E_c, nu=0.2) )
set_mat(dom.elems[:joints], MCJoint( E=E_c, nu=0.2, ft=f_t, mu=1.4, alfa=1, beta=1 ,wc=1.48e-4, ws=2.05e-5))  


# Acompanhamento de nó
noinf1 = dom.nodes[:(x==($Lviga/2-$lfiss) && y==0 && z==0)]
nodedat1 = NodeTracker(noinf1)
set_trackers(dom, nodedat1)

noinf2 = dom.nodes[:(x==($Lviga/2+$lfiss) && y==0 && z==0)]
nodedat2 = NodeTracker(noinf2)
set_trackers(dom, nodedat2)

#nosup = dom.nodes[:(x==$Lviga/2 && y==$hviga+$hapoio && z==$Eviga/2)][1]
facesup = FacesTracker( dom.faces[:(y==$hviga+$hapoio)] )
set_trackers(dom, facesup)

# Condicoes de contorno e solucao
bc1 = EdgeBC(:(x==($Lviga-$Lparcial)/2 && y==0), ux=0, uy=0, uz=0)
bc2 = EdgeBC(:(x==(($Lviga-$Lparcial)/2+$Lparcial) && y==0), uy=0)

#nocentro = dom.nodes[:(x==($Lviga/2) && y==$hviga+$hapoio && z==$Eviga/2)][1]
bc3 =FaceBC(dom.faces[:(y==$hviga+$hapoio)], uy=-w)

set_bc(dom, bc1, bc2, bc3)

solve!(dom, auto=true, nincs=10, scheme="ME", precision=1e-1, savesteps=true)
save(nodedat1, "noinf1.dat")
save(nodedat2, "noinf2.dat")
save(facesup, "facesup.dat")
