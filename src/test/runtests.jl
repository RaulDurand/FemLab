path  = dirname(@__FILE__)
tests = readdir(path)

for t in tests
    if t[1:2]=="t_"
        try
            println()
            println("Running test ", t, ":")
            include(t)
        catch
            warn("test", t, "failed.")
        end
    end
end
