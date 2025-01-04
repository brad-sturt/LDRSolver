using Distances
using LinearAlgebra
using Random
using SparseArrays
using JuMP
using Gurobi

# We consider the production problem from Ben-Tal et al (2004). 

include("../src/lp_solvers.jl")
include("../src/algorithm_full.jl")
include("../src/bental_setup.jl")

function main()

    ###########################################################################
    # Run the active set method for 
    # E=[10,20,30,40,50] and time_factor=10.
    ###########################################################################

    # Create output file
    outfile = string("../output/figure_2.csv")
    f_out = open(outfile,"w")

    # Create header to the output file
    row_string = string("E", ",",           # Number of factories
                        "time_factor", ",", # Time factor
                        "removal", ",",     # true if we wil remove at end of each iteration
                        "T", ",",           # Number of time periods
                        "iteration", ",",   # Iteration
                        "obj_val", ",",     # Current objective value
                        "A", ",",           # Length of current active set
                        "y_nonzero", ",",   # Number of nonzero parameters in y
                        "time", ",",        # Time for current iteration
                        "zeta_nonzero",",", # Number of nonzero parameters in ζ
                        "M", ",",           # big M
                        "hitting_M"         # true if there exists λ[i] = M
                        )
    row_string = string(row_string, "\n")
    print(f_out,row_string)
    flush(f_out)

    # Run experiment
    for E=[10,20,30,40,50], time_factor=[10],removal=[true]
       
        ###########################################################################
        # Set up the problem instance
        ###########################################################################
    
        a,b,c,D_min,D_max = GetProductionProblem(E,time_factor)


        ###########################################################################
        # Run algorithm and write to output
        ###########################################################################

        time_limit = 300
        ActiveSetMethod(a,b,c,D_min,D_max,f_out,removal,time_limit)
    end

    close(f_out)
end

main()
