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
    # Run the various methods for E=5 and time_factor=[2,4,6,8,10]. This includes:
    #
    #  - Active set method 
    #  - Barrier method (to 10%, 1%, and 0.1% accuracy)
    ###########################################################################

    for solution_method in ["dual", "activeset"]
        
    
        ############################################ 
        # Run the various algorithms
        ############################################ 

        if solution_method == "dual"
            
            println("=============== Dual ===============")

            # Create output file
            outfile = string("../output/figure_1_barrier.csv")
            f_out = open(outfile,"w")

            # Create header to the output file
            row_string = string("E", ",",           # Number of factories
                                "time_factor", ",", # Time factor
                                "T", ",",           # Number of time periods
                                "obj_val", ",",     # Current objective value
                                "y_nonzero", ",",   # Number of nonzero parameters in y
                                "time", ",",        # Time for current iteration
                                "zeta_nonzero", ",",
                                "Method", ",",
                                "gap") # Number of nonzero parameters in ζ)
            row_string = string(row_string, "\n")
            print(f_out,row_string)
            flush(f_out)

            # Run experiment
            for E=[5], time_factor=[2,4,6,8,10], method = [2], gap = [0.1,0.01,0.001]


                ###########################################################################
                # Set up the problem instance
                ###########################################################################
            
                a,b,c,D_min,D_max = GetProductionProblem(E,time_factor)


                ###########################################################################
                # Run algorithm
                ###########################################################################

                time_elapsed =  @elapsed  obj_val_D, y_D, ζ_D, λ_D = DualLPSolver(a,b,c,D_min,D_max,method,gap)


                ###########################################################################
                # Write to output
                ###########################################################################

                row_string = string(E, ",",           # Number of factories
                                    time_factor, ",", # Time factor
                                    time_factor*24, ",",           # Number of time periods
                                    obj_val_D, ",",     # Current objective value
                                    count(!iszero, y_D), ",", # Number of nonzero parameters in y
                                    time_elapsed, ",",        # Time for current iteration
                                    count(!iszero, Array(ζ_D)), ",",
                                    method,",",
                                    gap) 
                row_string = string(row_string, "\n")
                print(f_out,row_string)
                flush(f_out)

            end
            close(f_out)

        elseif solution_method == "activeset"
            println("=============== New Algorithm ===============")

            # Create output file
            outfile = string("../output/figure_1_active_set.csv")
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
            for E=[5], time_factor=[2,4,6,8,10],removal=[true]
               
                ###########################################################################
                # Set up the problem instance
                ###########################################################################
            
                a,b,c,D_min,D_max = GetProductionProblem(E,time_factor)


                ###########################################################################
                # Run algorithm and write to output
                ###########################################################################

                ActiveSetMethod(a,b,c,D_min,D_max,f_out,removal)
            end

            close(f_out)
        else
            error("Unknown method type, ", solution_method)
        end
    end
end

main()
