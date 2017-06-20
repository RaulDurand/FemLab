import Base.copy!

# Abstract type for material
# ==========================

"""
`Material`

Abstract type for objects used to store material parameters and to define
the bejaviour of elements.
"""
abstract type Material end

function copy!(target::Material, source::Material)
    @show("Trial function: This sould not be used")
    # source and target must be of the same type
    T = typeof(source)
    for fld in fieldnames(source)
        if filedtype(T, fld) <: Array
            if isdefined(target,fld)
                getfield(target, fld)[:] = getfield(source, fld)[:]
            else
                setfield!(target, fld, copy(getfield(source, fld)))
            end
        else
            setfield!(target, fld, getfield(source, fld))
        end
    end
end


# Reads material parameters from a json file
function read_prms(filename::String)

    # read file
    file = open(filename, "r")
    data = JSON.parse(file)
    close(file)

    # parse materials
    mats_prms = Dict{String, Any}()
    for d in data
        name = d["name"]
        keys = d["prms"]
        vals = d["vals"]
        prms = Dict{Symbol, Float64}( Symbol(k) => v for (k,v) in zip(keys, vals) )
        mats_prms[name] = prms
    end

    return mats_prms
end


