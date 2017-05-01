# Constants
using FemLab
using Base.Test

path  = dirname(@__FILE__)
tests = readdir(path)
verbose = false

print_with_color(:green, "\x1b[1m", "\nRunning tests...\n", "\x1b[0m")

@testset begin
    for t in tests
        if length(t)<5; continue end
        if t[1:5]!="test_"; continue end

        print_with_color(:white, "Running file ", t,"...\n")
        include(t)
        println()
    end
end
