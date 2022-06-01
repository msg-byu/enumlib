#!/usr/local/bin/julia
origlines = readlines(open("rankTest.in"))

k = 2:3 # Different numbers of elements to sweep over
labels = ["/"*string(i-1) for i ∈ k] # label list for input file
lattices = ["SimpleCubic","BCC","FCC","HCP"]
basis = [["1.0 0.0 1.0", "0.0 1.0 0.0", "0.0 0.0 1.0"],
         ["1.0 1.0 -1.","1.0 -1. 1.0","-1. 1.0 1.0"],
         ["1.0 0.0 1.0", "0.0 1.0 0.0", "0.0 0.0 1.0"],
         ["0.5 0.8660254037844386 0.0", 
         "-0.5 0.8660254037844386 0.0",
         "0.0 0.0 1.632993161855452"]]
for (j,lat) ∈ enumerate(lattices)
    lines = deepcopy(origlines)
    lines[3:5]= basis[j]
    for (i,s) in enumerate(k)
        lines[6] = string(s)
        lines[8] = "0 0 0 0"*join(labels[1:i])
        idx = lattices[j]*"_"*string(s)
        infile = "rankTest_"*idx*".in"
        open(infile,"w") do f
            for i ∈ lines
                println(f,i)
            end
        end
        mv(infile,"struct_enum.in",force=true)
        cmd = "../src/enum.x"
        run(`$cmd`)
        mkpath(idx)
        println(pwd())
        mv("struct_enum.out",idx*"/struct_enum.out",force=true)
        cd(idx)
        elements = "Pd "^(i+1)
        run(`python ../makeStr.py all -species $elements`)
        #run(cmd)
        cd("..")    
    end

end
