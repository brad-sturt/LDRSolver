using JuMP
using Gurobi
using LinearAlgebra
using TickTock
using TimerOutputs
include("utils.jl")
gurobi_env = Gurobi.Env()

# A problem is defined by the following parameters
# - a[i] is an T-by-n matrix
# - b[i] is an T-by-1 vecor
# - c is a m-by-1 vector

# The primal solution is defined by the following parameters:
# - y[t,s,j] for all 1 ≤ s ≤ t ≤ T and 1 ≤ j ≤ n

# The dual solution is defined by the following parameters:
# - λ[0],…,λ[m] ≥ 0
# - ζ[i,s] for all 0 ≤ i ≤ m and 1 ≤ s ≤ T


# The following function solves the problem as an LP in the dual formulation
function DualLPSolver(a,b,c,D_min,D_max,method,gap=0,time_limit=1e9)

    # Get the problem dimensions
    m = length(a)-1
    T,n = size(a[0])

    # Get the nonzero indices
    nonzeros_a, nonzeros_b, nonzeros_a_tj,  nonzeros_a_it = GetNonzeros(a,b)
    
    
    # Construct the optimization problem
    model = direct_model(Gurobi.Optimizer(gurobi_env))::JuMP.Model
    set_silent(model)

    # Uncomment if want to terminate method early
    set_time_limit_sec(model, time_limit)

    # Set gurobi settings
    if method != -1
        set_optimizer_attribute(model, "Method", method)
        if method == 2
            set_optimizer_attribute(model, "Crossover", 0)
            if gap != 0
                set_optimizer_attribute(model, "BarConvTol", gap)
            end
        end
    end

    # Construct decision variables
    @variable(model, λ[i=0:m] ≥ 0)
    @variable(model, ζ[i=0:m,s=1:T])

    # Construct constraints
    @constraint(model, D1[t=1:T,s=1:T,j=1:n ; s ≤ t], sum(nonzeros_a_tj[(t,j)][i]*ζ[i,s] for i in keys(nonzeros_a_tj[(t,j)])) == 0)
    @constraint(model, D2a[i=0:m,s=1:T], ζ[i,s] ≤ D_max[s]*λ[i])
    @constraint(model, D2b[i=0:m,s=1:T], ζ[i,s] ≥ D_min[s]*λ[i])
    @constraint(model, λ[0] == 1)

    # Construct objective
    @objective(model, Max, - sum(c[i]*λ[i] for i=1:m) - sum(nonzeros_b[(i,s)] * ζ[i,s] for (i,s) in keys(nonzeros_b)))

    # Solve optimization problem
    optimize!(model)

    # Parse output
    @show num_constraints(model; count_variable_in_set_constraints=false)
    @show num_variables(model)
    if termination_status(model) == OPTIMAL
        println("Solution is optimal")
    #elseif termination_status(model) == TIME_LIMIT && has_values(model)
    elseif has_values(model)
        println("Solution is suboptimal, but a primal solution is available")
    else
        println("The model was not solved correctly.")
        return Inf, [], [], []
    end
    if primal_status(model) == FEASIBLE_POINT
        println("Found feasible primal soln")
    end
    if dual_status(model) == FEASIBLE_POINT
        println("Found feasible dual soln")
    end

    # Return the optimal objective value and optimal solutions
    return objective_value(model), dual.(D1), value.(ζ), value.(λ)
end 


