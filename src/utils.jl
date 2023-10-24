


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



