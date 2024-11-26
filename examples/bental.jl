using Distances
using LinearAlgebra
using Random
using SparseArrays
using JuMP
using Gurobi

# We consider the production problem from Ben-Tal et al (2004). 

include("../src/lp_solvers.jl")
include("../src/algorithm_full.jl")



# E     -- Number of factories
# T     -- Number of time periods evenly spaced over one calendar year
# Dmin  -- Minmum demand in each period, where Dmin[t] is minimum demand
#          on period t
# Dmax  -- Maximum demand in each period, where Dmax[t] is maximum demand
#          on period t
# α     -- Production costs, where α[t][e] is production cost for
#          factory e on period t
# p     -- Maximum production per period per factory, where p[t][e] is the
#          maximum production for factory e on period t
# Q     -- Maximum production per year per factory, where Q[e] is the maximum
#          production for factory e
# Vmin  -- Minimum inventory at warehouse per period
# Vmax  -- Maximum inventory at warehouse per period
# v1    -- Initial inventory

function _GetProductionProblem(E,T,Dmin,Dmax,α,p,Q,Vmin,Vmax,v1)



    ###########################################################################
    # Specify the dimensions of the multi-stage robust linear optimization
    # problem.
    ###########################################################################

    # Because we have to make decisions in each period without knowledge of
    # the current demand, and because we have an offset term, there will be
    # a total of T+1 periods, in which decisions are made in periods 
    # 1,…,T and uncertainty is revealed on periods 2,…,T. The uncertainty
    # on the first period is the constant 1.

    # Store the number of constraints
    m = E + 2*E*T + 2*T

    # Store the number of decision variables in each period
    n = E


    ###########################################################################
    # Initialize the problem and uncertainty set
    ###########################################################################
   
    a = Dict{Int64,SparseMatrixCSC}(i => spzeros(Float64,T+1,n) for i=0:m)
    b = Dict{Int64,SparseVector}(i => spzeros(Float64,T+1) for i=0:m)
    c = zeros(m)
  
 
    ###########################################################################
    # Populate the problem constraints
    ###########################################################################
   
    # Initialize the starting row 
    i = 1

    # Add the constraint that we cannot transport more inventory out of each 
    # location than is availablek
    for e=1:E
        for t = 1:T
            a[i][t,e] = 1
        end
        c[i] = Q[e]
        i += 1
    end

    # Add the nonnegativity constraint for each factory on each period
    for e=1:E
        for t = 1:T
            a[i][t,e] = -1
            i += 1
        end
    end
 
    # Add the upper bound constraint for each factory on each period
    for e=1:E
        for t = 1:T
            a[i][t,e] = 1
            c[i] = p[t][e]
            i += 1
        end
    end
   
    # Add the upper bound constraint for warehouse in each period
    for t=1:T
        for s=1:t,e=1:E
            a[i][s,e] = 1
        end
        for s=2:(t+1)
            b[i][s] = 1
        end
        c[i] = Vmax - v1
        i += 1
    end
 
    # Add the lower bound constraint for warehouse in each period
    for t=1:T
        for s=1:t,e=1:E
            a[i][s,e] = -1
        end
        for s=2:(t+1)
            b[i][s] = -1
        end
        c[i] = -Vmin + v1
        i += 1
    end


    ###########################################################################
    # Populate the problem objective function
    ###########################################################################
   
    for t=1:T,e=1:E
        a[0][t,e] = α[t][e]
    end


    ###########################################################################
    # Return the problem and uncertainty sets
    ###########################################################################

    D_min = Vector{Float64}()
    D_max = Vector{Float64}()
    push!(D_min,1)
    push!(D_max,1)
    for t=1:T
        push!(D_min,Dmin[t])
        push!(D_max,Dmax[t])
    end
    return a,b,c,D_min,D_max
end


function GetProductionProblem(E,time_factor)
    T = 24*time_factor          # Number of time periods
    d_nom = 1000/time_factor*[1 + 0.5*sin(π*(t-1)/12)  for t in 1:T]  # Nominal demand
    θ = 0.20        # Uncertainty level
    Dmin = d_nom*(1-θ)
    Dmax = d_nom*(1+θ) 

    
    # Production costs
    α_nom = [1.0 + (e-1) / (E - 1) for e in 1:E]
    α = [[α_nom[e] * (1 + 0.5*sin(π*(t-1)/12))  for e in 1:E] for t in 1:T]
    
    # Production bounds
    p = [[567/time_factor for e in 1:E] for t=1:T]  # Maximimum production per period
    Q = [13600 for e in 1:E]       # Maximumum production over all periods
    
    # Bounds on warehouse
    Vmin = 500      # Minimum inventory at warehouse
    Vmax = 2000     # Maximum inventory at warehouse
    v1 = Vmin       # Initial inventory (not provided in paper)
    
    # Create the problem instance
    return _GetProductionProblem(E,T,Dmin,Dmax,α,p,Q,Vmin,Vmax,v1)
end

function main(methods=["dual","activeset"])

    ###########################################################################
    # Run the various methods
    ###########################################################################

    for method in methods
        
    
        ############################################ 
        # Run the various algorithms
        ############################################ 

        if method == "dual"
            
            println("=============== Dual ===============")

            # Create output file
            outfile = string("../output/dual.csv")
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

                # Get number of nonzero parameters in y
                num_nonzero_y_D = count(!iszero, y_D)
                num_nonzero_ζ_D = count(!iszero, Array(ζ_D))

                row_string = string(E, ",",           # Number of factories
                                    time_factor, ",", # Time factor
                                    time_factor*24, ",",           # Number of time periods
                                    obj_val_D, ",",     # Current objective value
                                    num_nonzero_y_D, ",", # Number of nonzero parameters in y
                                    time_elapsed, ",",        # Time for current iteration
                                    num_nonzero_ζ_D, ",",
                                    method,",",
                                    gap) 
                row_string = string(row_string, "\n")
                print(f_out,row_string)
                flush(f_out)

            end
            close(f_out)


        elseif method == "activeset"
            println("=============== Active Set Method ===============")

            # Create output file
            outfile = string("../output/activeset.csv")
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
            error("Unknown method type, ", method)
        end
    end
end

