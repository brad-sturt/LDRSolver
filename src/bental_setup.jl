using Distances
using LinearAlgebra
using Random
using SparseArrays


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

