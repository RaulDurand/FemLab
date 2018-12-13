
export Beam

mutable struct BeamIpState<:IpState
    ndim::Int
    #σ::Float64
    #ε::Float64
    function BeamIpState(ndim::Int=3)
        return new(ndim)
    end
end

mutable struct Beam<:AbsBeam
    E::Float64
    A::Float64
    I::Float64
   xi::Float64
   ro::Float64

    function Beam(prms::Dict{Symbol,Float64})
        return  Beam(;prms...)
    end

    function Beam(;E=NaN, A=NaN, I=NaN, xi=NaN, ro=NaN, Ix=NaN, Iy=NaN, Iz=NaN)
        @assert E>0.0
        @assert A>0.0
        @assert I>0.0
        @assert xi>0.0
        @assert ro>0.0
        this = new(E,A,I,xi,ro)
        return this
    end
end

# Create a new instance of Ip data
new_ip_state(mat::Beam, ndim::Int) = BeamIpState(ndim)

#function set_state(ipd::BeamIpState, σ=NaN, ε=NaN)
    #if !isnan(σ); ipd.σ = σ end
    #if !isnan(ε); ipd.ε = ε end
#end

function calcT(mat::Beam, elem::Element, C)
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,1])/L
    return 

end

function elem_stiffness(mat::Beam, elem::Element)
    C  = elem_coords(elem)
    L  = norm(C[2,:]-C[1,:])
    L2 = L*L
    L3 = L*L*L
    EA = mat.E*mat.A
    EI = mat.E*mat.I

    K0 = [ EA/L     0         -EA/L      0      0          0
           0       12*EI/L3   6*EI/L2    0     -12*EI/L3   6*EI/L2
           0        6*EI/L2   4*EI/L     0      -6*EI/L2   2*EI/L 
          -EA/L     0          0         EA/L     0        0
           0      -12*EI/L3  -6*EI/L2    0      12*EI/L3  -6*EI/L2
           0        6*EI/L2   2*EI/L     0      -6*EI/L2   4*EI/L  ]


    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,1])/L
    T = eye(6)*c
    T[3,3] = T[6,6] = 1.0
    T[1,1] = T[2,2] = T[4,4] = T[5,5] = c
    T[1,2] = T[4,5] =  s
    T[2,1] = T[5,4] = -s
    return T'*K0*T
end           

function elem_mass(mat::Beam, elem::Element)
    C  = elem_coords(elem)
    L  = norm(C[2,:]-C[1,:])
    L2 = L*L
    L3 = L*L*L
    EA = mat.E*mat.A
    EI = mat.E*mat.I
    roAL = mat.ro*mat.A*L
    
    K0 = (roAL/420) * [ 140       0          0       70         0           0
                          0     156       22*L        0        54       -13*L
                          0     22*L      4*L2        0      13*L       -3*L2 
                         70       0          0      140         0           0
                          0      54       13*L        0       156       -22*L
                          0    -13*L     -3*L2        0     -22*L        4*L2]


    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,1])/L
    T = eye(6)*c
    T[3,3] = T[6,6] = 1.0
    T[1,1] = T[2,2] = T[4,4] = T[5,5] = c
    T[1,2] = T[4,5] =  s
    T[2,1] = T[5,4] = -s
    return T'*K0*T
end   

function ip_state_vals(mat::Beam, ipd::BeamIpState)
    return Dict()
end
