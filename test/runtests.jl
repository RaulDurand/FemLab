# Constants
using FemLab
using FactCheck

path  = dirname(@__FILE__)
tests = readdir(path)

#println(GREEN, BOLD, "\nRunning tests...", DEFAULT)
print_with_color(:green, "\x1b[1m", "\nRunning tests...\n", "\x1b[0m")

for t in tests
    if t[1:5]=="test_"
        include(t)
    end
end

FactCheck.exitstatus()
