


###########################################################################
# Create data structures which store the nonzero entries of the matrices a
# and b
# 
#    nonzeros_a = {(i,t,j) => a[i][t][j] : a[i][t][j] is nonzero}
#    nonzeros_b = {(s,i) => b[i][s] : b[i][s] is nonzero}
#    nonzeros_a_tj[(t,j)] = {i => a[i][t][j] : a[i][t][j] is nonzero}
# 
###########################################################################

function GetNonzeros(a,b)

    # Get the problem dimensions
    m = length(a)-1
    T,n = size(a[0])

    # Initialize the data structures
    nonzeros_a = Dict{Tuple{Int64,Int64,Int64}, Float64}()
    nonzeros_b = Dict{Tuple{Int64,Int64}, Float64}()
    nonzeros_b_s = Dict{Int64, Dict{Int64,Float64}}()
    nonzeros_a_tj = Dict{Tuple{Int64,Int64},Dict{Int64,Float64}}()
    nonzeros_a_it = Dict{Tuple{Int64,Int64},Dict{Int64,Float64}}()
    for t=1:T,j=1:n
        nonzeros_a_tj[(t,j)] = Dict{Int64,Float64}()
    end
    for t=1:T,i=0:m
        nonzeros_a_it[(i,t)] = Dict{Int64,Float64}()
    end
    for s=1:T
        nonzeros_b_s[s] = Dict{Int64,Float64}()
    end


    # Populate the data structures
    for i=0:m
        rows,cols,vals = findnz(a[i])

        for index = 1:length(rows)
            t = rows[index]
            j = cols[index]
            nonzeros_a_tj[(t,j)][i] = vals[index]
            nonzeros_a_it[(i,t)][j] = vals[index]
            nonzeros_a[(i,t,j)] = vals[index]
        end

    end
    for i=0:m
        rows,vals = findnz(b[i])
        for index = 1:length(rows)
            s = rows[index]
            nonzeros_b[(i,s)] = vals[index]
            nonzeros_b_s[s][i] = vals[index]
        end
    end

    return nonzeros_a, nonzeros_b, nonzeros_a_tj, nonzeros_a_it, nonzeros_b_s
end


function ComputeIandS_new(A,a,b)


    ###########################################################################
    # Get the problem dimensions
    ###########################################################################

    m = length(a)-1
    T = length(a[0])
    n = length(a[0][1])

 
    ###########################################################################
    # Construct the nonzero data structures
    ###########################################################################

    nonzeros_a, nonzeros_b, nonzeros_a_tj, nonzeros_a_it, nonzeros_b_s = GetNonzeros(a,b)


    ###########################################################################
    # Initialize the data structures
    ###########################################################################

    I_A = Set{Int64}()
    push!(I_A,0)
    S_A = Dict{Int64,Set{Int64}}()
    S_A[0] = Set{Int64}()
    S_A_min = Dict{Int64,Set{Int64}}()
    S_A_min[0] = Set{Int64}()
    S_A_max = Dict{Int64,Set{Int64}}()
    S_A_max[0] = Set{Int64}()


    ###########################################################################
    # Populate the data structures
    ###########################################################################
 
    # Populate the entries for S_A
    for (t,s,j) in A
        for i in keys(nonzeros_a_tj[(t,j)])
            if i ∉ I_A
                push!(I_A,i)
                S_A[i] = Set{Int64}()
            end
            push!(S_A[i],s)
        end
    end
   
    # Add the nonzero entries for S_A_min and S_A_max
    for (i,s) in keys(nonzeros_b)
        if nonzeros_b[(i,s)] < 0
            if i ∉ I_A
                push!(I_A,i)
                S_A_min[i] = Set{Int64}()
            end
            push!(S_A_min[i],s)
        else
            if i ∉ I_A
                push!(I_A,i)
                S_A_max[i] = Set{Int64}()
            end
            push!(S_A_max[i],s)
        end
    end



    ###########################################################################
    # Return the data structures
    ###########################################################################
   
    return I_A, S_A, S_A_min, S_A_max

end


function ComputeIandS(A,a,b)


    ###########################################################################
    # Get the problem dimensions
    ###########################################################################

    m = length(a)-1
    T = length(a[0])
    n = length(a[0][1])

 
    ###########################################################################
    # Construct the nonzero data structures
    ###########################################################################

    nonzeros_a, nonzeros_b, nonzeros_a_tj, nonzeros_a_it, nonzeros_b_s = GetNonzeros(a,b)


    ###########################################################################
    # Initialize the data structures
    ###########################################################################

    I_A = Set{Int64}()
    push!(I_A,0)
    S_A = Dict{Int64,Set{Int64}}()
    S_A[0] = Set{Int64}()


    ###########################################################################
    # Populate the data structures
    ###########################################################################
 
    # Populate the entries for S_A
    for (t,s,j) in A
        for i in keys(nonzeros_a_tj[(t,j)])
            if i ∉ I_A
                push!(I_A,i)
                S_A[i] = Set{Int64}()
            end
            push!(S_A[i],s)
        end
    end
    for (i,s) in keys(nonzeros_b)
        if i ∉ I_A
            push!(I_A,i)
            S_A[i] = Set{Int64}()
        end
        push!(S_A[i],s)
    end



    ###########################################################################
    # Return the data structures
    ###########################################################################
   
    return I_A, S_A

end

