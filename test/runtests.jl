# Constants
using FemLab
using FactCheck

path  = dirname(@__FILE__)
tests = readdir(path)
verbose = false

#println(GREEN, BOLD, "\nRunning tests...", DEFAULT)
print_with_color(:green, "\x1b[1m", "\nRunning tests...\n", "\x1b[0m")

for t in tests
    if length(t)<5; continue end
    if t[1:5]!="test_"; continue end

    print("Running file ", t,"...")
    include(t)
    println()
end

FactCheck.exitstatus()
