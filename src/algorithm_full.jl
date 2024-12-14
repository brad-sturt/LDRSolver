using JuMP
using Gurobi
using LinearAlgebra
using TickTock
using TimerOutputs
using RowEchelon

include("utils.jl")
gurobi_env = Gurobi.Env()


function ActiveSetMethod(a,b,c,D_min,D_max,f_out,removal=false,time_limit=1e9)

    
    ###########################################################################
    # Get the problem dimensions
    ###########################################################################
    
    m = length(a)-1
    T,n = size(a[0])


    ###########################################################################
    # Set the initial value of the big M in the uncertainty sets
    ###########################################################################
    
    # TODO Improve how we select this value
    M = 1e7
    counter = 0
    cum_time = 0


    ###########################################################################
    # Construct the nonzero data structures
    ###########################################################################

    nonzeros_a, nonzeros_b, nonzeros_a_tj, nonzeros_a_it, nonzeros_b_s = GetNonzeros(a,b)


    ###########################################################################
    # Create an initial active set A
    ###########################################################################
    
    A = Set{Tuple{Int64,Int64,Int64}}()


    ###########################################################################
    # Create a dictionary for keeping track of the objective value when we
    # added the tuple
    ###########################################################################

    A_info = Dict{Tuple{Int64,Int64,Int64}, Float64}()


    ###########################################################################
    # Initialize active set
    ###########################################################################

    # We initially populate A with all of the offsets
    for t=1:T, j=1:n
        push!(A,(t,1,j))
        A_info[(t,1,j)] = Inf
    end

    # We initially populate A with all of markovian rules
    for t=1:T, j=1:n
        push!(A,(t,t,j))
        A_info[(t,t,j)] = Inf
    end


    ###########################################################################
    # Create an array to store the random permutations, which are generated
    # in order to select which tuples to add into active set in each iteration
    #
    # NOTE: These arrays were used when we tried using randomization to
    #       decide which tuples should be added into the active set
    ###########################################################################

    random_periods = zeros(Int,T)


    ###########################################################################
    # Create a temporary variable for the return variables
    ###########################################################################
    
    obj_val = Inf
    obj_val_old = Inf
    y_opt = 0
    num_iterations_without_improvement = 0


    ###########################################################################
    # Create a sparse matrix version of A, which will be useful in each
    # iteration when computing the violated constraints
    # 
    # NOTE: A_sparse[t] is an (m+1)x(n) sparse matrix
    ###########################################################################

    a_sparse = [spzeros(Float64,m+1,n) for t=1:T]
    for (i,t,j) in keys(nonzeros_a)
        a_sparse[t][i+1,j] = nonzeros_a[i,t,j]
    end


    ###########################################################################
    # Perform iterations
    ###########################################################################

    while true
    
        ###########################################################################
        # Update counter
        ###########################################################################

        counter += 1
    
        temp_time = @elapsed λ_opt,y_opt, obj_val, ζ_sparse = _SolveRestrictedProblem(
                                    a,b,c,D_min,D_max,A,M)


        ###################################################################
        # Determine which tuples to add to active set
        ###################################################################

        # Create a variable to keep track of whether we added a tuple
        added_tuple = false

        # Create the temp matrix
        temp =  [dropzeros(a_sparse[t]'*ζ_sparse[s]) for t=1:T,s=1:T]
        nonzeros_temp = Dict{Tuple{Int64,Int64,Int64},Float64}()
        for t=1:T, s=1:t
            rows,vals = findnz(temp[t,s])
            for index = 1:length(rows)
                j = rows[index]
                if abs(vals[index]) > 1e-4
                    nonzeros_temp[(t,s,j)] = vals[index]
                end
            end
        end
        
        # Iterate over decision variables
        random_decisions = randperm(n)
        for j in random_decisions, t in 1:T
            
            # We will now choose the tuples to add to I for 
            # each decision variable x_{t,j}(⋅).
            # Currently, we will just add all of the tuples whose
            # corresponding constraints are violated (that is,
            # temp[t,s][j,ℓ] != 0).
            # Create a random ordering of past periods and
            # uncertain variables in each past period
            # (The randomness is currently unnecessary)
            randperm!(random_periods)

            # Iterate over the random periods and uncertanties
            for s=random_periods
                # If the selection is nonzero, skip
                if (t,s,j) ∉ keys(nonzeros_temp)
                    continue
                end
                # Check if the coefficient should enter
                if abs(nonzeros_temp[t,s,j]) > 1e-4
                    
                    # Add the tuple 
                    added_tuple = true
                    push!(A, (t,s,j))
                    A_info[(t,s,j)] = obj_val
                    break
                end
            end
        end


        ###################################################################
        # Determine whether we have a feasible primal linear decision rule
        ###################################################################

        hitting_M = false
        for i=1:m
            if λ_opt[i] ≥ M - 1e-3
                hitting_M = true
                break
            end
        end


        ###################################################################
        # Determine which tuples to remove from A
        ###################################################################

        if removal && !hitting_M

            # Iterate over each of the tuples in I
            for (t,s,j) in A

                # See if the tuple is a candidate for removal by
                # by checking if Y_opt[t,s][j,ℓ] is equal to zero
                if abs(y_opt[t][s,j]) ≤ 1e-8

                    # Check if the objective value has changed since
                    # the tuple (t,s,j,ℓ) was added to I
                    if obj_val < A_info[(t,s,j)] - 1e-8

                        # Remove the tuple from active set
                        delete!(A, (t,s,j))
                        delete!(A_info, (t,s,j))
                    end
                end
            end
        end


        ###################################################################
        # Give useful debugging information
        ###################################################################

        cum_time += temp_time

        println("Parameters: ", sum(sum(sum(1 for j=1:n) for s=1:t) for t=1:T), 
                ";\tLength of A: ", length(A), 
                ";\t ||y||0: ", sum(count(!iszero, y_opt[t]) for t=1:T),
                ";\t Objective value: ",  obj_val,
                ";\t Solver time: ", temp_time)


        ###################################################################
        # Log to file
        ###################################################################

        row_string = string(n, ",",                                       # Number of factories
                            Int((T-1)/24), ",",                           # Time factor
                            removal, ",",                                 # true if we will remove at end of each iteration
                            T-1, ",",                                     # Number of time periods
                            counter, ",",                                 # Iteration
                            obj_val, ",",                                 # Current objective value
                            length(A), ",",                               # Length of current active set
                            sum(count(!iszero, y_opt[t]) for t=1:T), ",", # Number of nonzero parameters in y
                            temp_time, ",",                               # Time for current iteration
                            sum(count(!iszero, 
                                      ζ_sparse[t]) for t=1:T), ",",       # Number of nonzero parameters in ζ
                            M, ",",                                       # big M
                            hitting_M                                     # true if there exists λ[i] = M
                            )
        row_string = string(row_string, "\n")
        print(f_out,row_string)
        flush(f_out)


        ###################################################################
        # Check if new objective value is different
        ###################################################################
        
        if obj_val < obj_val_old - 1e-3
            num_iterations_without_improvement = 0
        else
            num_iterations_without_improvement += 1
        end
        obj_val_old = obj_val


        ###################################################################
        # Break if we didn't update I. Otherwise, do another iteration
        ###################################################################

        if !added_tuple || cum_time > time_limit || num_iterations_without_improvement ≥ 100
            break
        else
            continue
        end
    end
   
    println("-----------------------------------------------------------------")
end

function _SolveRestrictedProblem(a,b,c,D_min,D_max,A,M)
    
    ###########################################################################
    # Get the problem dimensions
    ###########################################################################
    
    m::Int64 = length(a)-1
    T::Int64,n::Int64 = size(a[0])


    ###########################################################################
    # Construct the nonzero data structures
    ###########################################################################

    nonzeros_a, nonzeros_b, nonzeros_a_tj, nonzeros_a_it, nonzeros_b_s = GetNonzeros(a,b)


    ###########################################################################
    # Initialize the restricted optimization problem
    ###########################################################################

    model = direct_model(Gurobi.Optimizer(gurobi_env))::JuMP.Model
    set_silent(model)


    ###########################################################################
    # Compute the set of unique tuples for the given active set
    ###########################################################################

    # Let A_s_to_tj[s] be the tuples (t,j) for which (t,s,j) ∈ A
    A_s_to_tj = Dict{Int64,Set{Tuple{Int64,Int64}}}()
    for s=1:T
        A_s_to_tj[s] = Set{Tuple{Int64,Int64}}()
    end
    for (t,s,j) ∈ A
        push!(A_s_to_tj[s],(t,j))      
    end

    # Let I_A be the nT-by-nT^2 matrix
    I_A_rows = Vector{Int64}()
    I_A_cols = Vector{Int64}()
    I_A_vals = Vector{Float64}()
    for s=1:T
        for (t,j) in A_s_to_tj[s]
            push!(I_A_rows,(t-1)*n+j)
            push!(I_A_cols,(s-1)*T*n+(t-1)*n+j)
            push!(I_A_vals, 1)
        end
    end
    I_A = sparse(I_A_rows,I_A_cols,I_A_vals,T*n,T^2*n)

    # Let D be the matrix
    D_rows = Vector{Int64}()
    D_cols = Vector{Int64}()
    D_vals = Vector{Float64}()
    for (i,t,j) in keys(nonzeros_a)
        push!(D_rows, i+1)
        push!(D_cols, (t-1)*n+j)
        push!(D_vals, nonzeros_a[i,t,j])
    end 
    D = sparse(D_rows,D_cols,D_vals,m+1,T*n)

    # Do the multiplication
    temp = D*I_A
    rows,cols,vals = findnz(temp)

    # Get the unique values
    unique_vals = Set{Float64}()
    for v in vals
        push!(unique_vals,v)
    end
    for v in values(nonzeros_b)
        push!(unique_vals,v)
    end


    # Create dictionary for each unique value

    # map[v][(t-1)*n+j,k] is the bucket for rows
    # whose value in column T*(s-1)*n + (t-1)*n + j is
    # v and whose previous bucket was k. Note that 
    # s = bucket_updated[v][(t-1)*n+j,k]
    map = Dict{Float64,Dict{Tuple{Int64,Int64},Int64}}(v => Dict{Tuple{Int64,Int64},Int64}() for v in unique_vals)

    # Create dictionary for when each unique value was updated
    map_updated = Dict{Float64,Dict{Tuple{Int64,Int64},Int64}}(v => Dict{Tuple{Int64,Int64},Int64}() for v in unique_vals)

    # Create buckets dictionary
    bucket = Dict{Tuple{Int64,Int64},Int64}()
    for i=0:m,s=1:T
        bucket[i,s] = 1
    end

    # Populate the buckets
    num_buckets = 1 
    index = 1
    s_iteration = 1
    s::Int64 = 0
    v = 0
    i::Int64 = 0
    t::Int64 = 0
    temp_index::Int64 = 0
    j::Int64 = 0
    k::Int64 = 0

    while true
        if s_iteration > T
            break
        end

        if index > length(rows) || ((cols[index]-1) ÷ (T*n)) + 1 > s_iteration

            # Complete the iteration by accounting for the b vector
            for i in keys(nonzeros_b_s[s_iteration])
                k = bucket[i,s_iteration]
                v = nonzeros_b_s[s_iteration][i]
                if (0,k) ∉ keys(map[v])
                    num_buckets += 1
                    map[v][0,k] = num_buckets
                    map_updated[v][0,k] = s_iteration
                elseif map_updated[v][0,k] != s_iteration
                    num_buckets += 1
                    map[v][0,k] = num_buckets
                    map_updated[v][0,k] = s_iteration
                end
                bucket[i,s_iteration] = map[v][0,k]
            end

            # Reset num_buckets
            num_buckets = 1
            
            # Skip to the next s_iteration
            s_iteration += 1

            # Go back to the beginning
            continue
        end

        s = ((cols[index]-1) ÷ (T*n)) + 1
        v = vals[index]
        i = rows[index] - 1
        t = ((cols[index] - (s-1)*T*n - 1) ÷ n)+1
        j = cols[index] - (s-1)*T*n - (t-1)*n
        k = bucket[i,s]
        temp_index = (t-1)*n+j
        if (temp_index,k) ∉ keys(map[v])
            num_buckets += 1
            map[v][temp_index,k] = num_buckets
            map_updated[v][temp_index,k] = s
        elseif map_updated[v][temp_index,k] != s
            num_buckets += 1
            map[v][temp_index,k] = num_buckets
            map_updated[v][temp_index,k] = s
        end
        bucket[i,s] = map[v][temp_index,k]

        index += 1
    end

    # Perform a reindexing of the buckets so that the first
    # bucket starts at index 1, the second bucket starts at 
    # index 2, and so on.
    # unique_buckets[s] = unique bucket indices for period s
    unique_buckets = Dict{Int64,Set{Int64}}()
    for s=1:T
        unique_buckets[s] = Set{Int64}()
        for i=0:m
            push!(unique_buckets[s],bucket[i,s])
        end
    end
    # reindex[s,k] = new bucket index for k
    K = Vector{Int64}()
    for s=1:T
        push!(K,  0)
    end

    reindex = Dict{Tuple{Int64,Int64},Int64}()
    for s=1:T
        for k in unique_buckets[s]
            K[s] += 1
            reindex[s,k] = K[s]
        end
    end

    new_bucket = Dict{Tuple{Int64,Int64},Int64}()
    KS_A = Set{Tuple{Int64,Int64}}()
    for i=0:m
        for s=1:T
            new_bucket[i,s] = reindex[s,bucket[i,s]]

            # If bucket[i,s] == 1, then there are only nonzeros
            # in the row, and so we do not include it in KS_A
            if bucket[i,s] != 1
                push!(KS_A, (new_bucket[i,s],s))
            end
        end
    end

    # Replace bucket with new_bucket
    bucket = new_bucket

    # Find a representative row for each k
    representative_row = Dict{Tuple{Int64,Int64},Int64}()
    for s=1:T
        for i=0:m
            k = bucket[i,s]
            if (k,s) ∉ keys(representative_row)
                representative_row[k,s] = i
            end
        end
    end



    
    ###########################################################################
    # Create the decision variables
    ###########################################################################

    # Create the decision variables
    @variable(model, 0 ≤ λ[i=0:m] ≤ M)  
    max_K = maximum([K[s] for s=1:T])
    @variable(model, ζ[k=1:max_K,s=1:T; k ≤ K[s]])

    
    ###########################################################################
    # Create the constraints
    ###########################################################################

  
    # Create a variable for keeping track of the number of constraints.
    # This variable will be used for efficiently extracting the optimal
    # dual solution. 
    num_constraints = 0

    # Create a map to store a mapping from the constraint to the
    # tuple (t,s,j)
    constraints_map = Dict{Int64,Tuple{Int64,Int64,Int64}}()
    constraints_offset_map = Dict{Int64,Tuple{Int64,Int64}}()
    constraints_vec = collect(A)
    
    # Add the equality constraints

    # Construct the set KS_A, where (k,s) ∈ KS_A if there
    # is a nonzero in a[representative_row[k,s]]

    # Determine the buckets k corresponding to each
    # of the (t,s,j)
    representative_row_is_to_k = Dict{Tuple{Int64,Int64},Int64}()
    for s=1:T
        for k=1:K[s]
            i = representative_row[k,s]
            representative_row_is_to_k[i,s] = k
        end
    end

    relevant_buckets = Dict{Tuple{Int64,Int64,Int64},Set{Int64}}()

    debug_counter = 0
    for (t,s,j) ∈ A
        relevant_buckets[t,s,j] = Set{Int64}()
        v1 = K[s]
        v2 = length(nonzeros_a_tj[t,j])
        if v1 < v2
            for k=1:K[s]
                debug_counter += 1
                i = representative_row[k,s]
                if i ∈ keys(nonzeros_a_tj[t,j])
                    push!(relevant_buckets[t,s,j],k)
                end
            end
        else
            for i in keys(nonzeros_a_tj[t,j])
                debug_counter += 1
                if (i,s) in keys(representative_row_is_to_k)
                    k = representative_row_is_to_k[i,s]
                    push!(relevant_buckets[t,s,j],k)
                end
            end
        end
    end


    # Create the initial constraint
    num_constraints += 1
    constraints_map[1] = constraints_vec[1]
    (t,s,j) = constraints_map[1]
    @constraint(model,
        EqualityConstraints[num_constraints:num_constraints],
        sum(ζ[k,s]*nonzeros_a_tj[t,j][representative_row[k,s]] for k ∈ relevant_buckets[t,s,j] ) == 0,
        container=SparseAxisArray)

    # Now add the remaining constraints
    for index=2:length(constraints_vec)
        num_constraints += 1
        constraints_map[index] = constraints_vec[index]
        (t,s,j) = constraints_map[index]
        EqualityConstraints[num_constraints] = @constraint(model,
        sum(ζ[k,s]*nonzeros_a_tj[t,j][representative_row[k,s]] for k ∈ relevant_buckets[t,s,j] ) == 0)
    end



    # Add the inequality constraints for ζ
    # Find the constraints that correspond to bucket k for each
    # period s
    bucket_to_constraints = Dict{Tuple{Int64,Int64},Set{Int64}}()
    for s=1:T
        for k=1:K[s]
            bucket_to_constraints[k,s] = Set{Int64}()
        end
    end
    for s=1:T
        for i=0:m
            push!(bucket_to_constraints[bucket[i,s],s], i)
        end
    end
    @constraint(model,
        InequalityConstraints_Max[s=1:T,k=1:K[s]; (k,s) ∈ KS_A],
        ζ[k,s] ≤ D_max[s]*sum(λ[i] for i ∈ bucket_to_constraints[k,s]),
        container=SparseAxisArray)
    @constraint(model,
        InequalityConstraints_Min[s=1:T,k=1:K[s]; (k,s) ∈ KS_A],
        -ζ[k,s] ≤ -D_min[s]*sum(λ[i] for i ∈ bucket_to_constraints[k,s]),
        container=SparseAxisArray)

    for s=1:T, k=1:K[s]
            num_constraints += 2
    end

    # NEW Add variable for λ_sum
    @variable(model, λ_sum)
    @constraint(model, λ_sum == sum(λ[i] for i=0:m))
    num_constraints += 1

    # Add the tuples that were missing from KS_A
    @constraint(model,
        InequalityConstraints_Max_Zeros[s=1:T,k=1:K[s]; (k,s) ∉ KS_A],
        ζ[k,s] ≤ D_max[s]*(λ_sum - sum(sum(λ[i] for i ∈ bucket_to_constraints[k_temp,s]) for k_temp=1:K[s] if k_temp != k)),
        container=SparseAxisArray)
    @constraint(model,
        InequalityConstraints_Zeros[s=1:T,k=1:K[s]; (k,s) ∉ KS_A],
        -ζ[k,s] ≤ -D_min[s]*(λ_sum - sum(sum(λ[i] for i ∈ bucket_to_constraints[k_temp,s]) for k_temp=1:K[s] if k_temp != k)),
        container=SparseAxisArray)


    # Add the constraint that λ[0] is equal to one
    @constraint(model, LambdaZero, λ[0] == 1)
    num_constraints += 1



    ###########################################################################
    # Create the objective
    ###########################################################################

    @objective(model,
        Max, 
        - sum( c[i]*λ[i] for i=1:m)
        - sum(sum( b[representative_row[k,s]][s]*ζ[k,s] for k=1:K[s] if b[representative_row[k,s]][s] != 0) for s=1:T))


    ###########################################################################
    # Solve restricted optimization problem
    ###########################################################################

    # Set Gurobi to use dual simplex
    set_optimizer_attribute(model, "Method", 1)

    # Solve the optimization problem
    optimize!(model)


    ###########################################################################
    # Extract information about the solver
    # NOTE: This information is not currently used
    ###########################################################################

    num_constraints_moi = MOI.get(model, Gurobi.ModelAttribute("NumConstrs"))
    num_variables_moi = MOI.get(model, Gurobi.ModelAttribute("NumVars"))
    num_simplex_iterations = MOI.get(model, Gurobi.ModelAttribute("IterCount"))

    
    ###########################################################################
    # Extract the optimal primal solution
    ###########################################################################


    # Initialize the data structures for the primal solution 
    ζ_sparse = [spzeros(Float64,m+1) for s=1:T]
    λ_opt_orig = value.(λ)
    ζ_opt_orig = value.(ζ)
    λ_opt = Dict{Int64,Float64}()
    ζ_opt = Dict{Tuple{Int64,Int64},Float64}()
    for i=0:m
        λ_opt[i] = λ_opt_orig[i]
    end
    for s=1:T
        for k=1:K[s]
            ζ_opt[k,s] = ζ_opt_orig[k,s]
        end
    end

    # Determine the denominators
    denominators = Dict{Int64,Vector{Float64}}()

    for s=1:T
        denominators[s] = zeros(K[s])
        for i=0:m
            k = bucket[i,s]
            denominators[s][k] += λ_opt[i]
        end
    end
    for s=1:T
        for i=0:m
            if λ_opt[i] !=0
                ζ_sparse[s][i+1] = ζ_opt[bucket[i,s],s]*λ_opt[i] / denominators[s][bucket[i,s]]
            end
        end
    end



    ###########################################################################
    # Extract the optimal dual solution
    ###########################################################################


    # Create a data structure to store the dual solution
    dual_solution = zeros(num_constraints)  

    # Extract the optimal dual solution 
    GRBgetdblattrarray(
        backend(model),
        "Pi",
        0,
        num_constraints,
        dual_solution)

    # Initialize the data structure for the linear decision rule
    y = [spzeros(Float64,T,n) for t=1:T]

    # Populate the data structure for the linear decision rule
    for k in keys(constraints_map)
        y[constraints_map[k][1]][constraints_map[k][2],constraints_map[k][3]] = dual_solution[k]
    end


    ###########################################################################
    # Return the optimal solution
    ###########################################################################

    return λ_opt, y, objective_value(model), ζ_sparse
end


