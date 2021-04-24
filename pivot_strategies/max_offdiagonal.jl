using LinearAlgebra

"""
Utility function for bparlett and rook pivoting strategy.
Search the max element in the subdiagonal and non-diagonal part of A matrix starting from column_index .
"""
function max_offdiagonal(A::AbstractMatrix, column_index :: Number)
    
    # n is the number of rows of matrix A
    n = size(A,1)
    max_value = 0
    
    # Return nothing if no subdiagonal values possible
    if (column_index >= n)
        max_value = 0
        max_value_index = 0
    else
        for i in 1:n     
           
            if !(i==column_index)
                v = abs(A[i, column_index]) 
                if(v > max_value) 
                    max_value = v
                    max_value_index = i
                end
            end

        end   
    end

    return max_value, max_value_index
end