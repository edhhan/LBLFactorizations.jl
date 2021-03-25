using LinearAlgebra

"""
Utility function for bparlett and rook pivoting strategy
"""
function max_subdiagonal(A::AbstractMatrix, column_index :: Number)
    #n is the number of rows of matrix A
    n = size(A,1)
    max_value = 0
    #Return nothing if no subdiagonal values possible
    if (column_index >= n)
        max_value = 0
        max_value_index = 0
    else
        for i in column_index+1:n
            #Present value
            v = abs(A[i,column_index])
            
            if(v > max_value) 
                max_value = v
                max_value_index = i
            end

        end   
    end

    return max_value, max_value_index

end