# Constants
using FemLab
using FactCheck

path  = dirname(@__FILE__)
tests = readdir(path)

println(GREEN, BOLD, "\nRunning tests...", DEFAULT)
for t in tests
    if t[1:2]=="t_"
        include(t)
    end
end

FactCheck.exitstatus()
