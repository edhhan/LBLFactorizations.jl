using LinearAlgebra
include("max_subdiagonal.jl")

"""
Partial pivoting strategy in O(n). Provides the pivot according to the bunch-kaufmann pivoting strategy.
A search is done in at most two columns (always the first one).
    
    The chapter 11 in [1] describe the following elements in the strategy:
        -ω_1 : the maximimum magnitude of any subdiagonal entry in column 1
        -r   : the row index of ω_1
        -ω_r : the maximimum magnitude of any off-diagonal entry in column r
        -α : a parameter derived by minimizing the bound on the element growth [1]

    Input :
        -A : a matrix 
    Output :
        -pivot_size: size of the pivot, i.e size ∈ {0, 1, 2}
        -pivot : the indices linked to the permutation (pivoting) 
                -if size==1 then pivot = tuple(idx1,idx2)
                -if size==2 then pivot = [tuple(idx1,idx2), tuple(idx3,idx4)] 
                -if size==0 then pivot = -1 (special case where we can factorized more efficiently a column for L)
     
Based on :     
[1] N. J. Higham, Accuracy and stability of numerical algorithms (Chapter 11). SIAM, 2002.
"""
function bkaufmann(A::AbstractMatrix{T}) where T
    α = (1+sqrt(17))/8  # from [1]
    n = size(A,1)
    pivot = []
    pivot_size = 0

    # FIRST SEARCH : looks for the largest absolute value in the subdiagonal element of column 1 
    ω_1, r = max_subdiagonal(A, 1)
    a_11 = abs(A[1,1])
    
    # Early stopping : no pivoting required and trivial solution for the factorization (special case)
    if (ω_1 == 0)
        pivot_size = 0
        push!(pivot, -1)
        return pivot, pivot_size
    end

    # The pivot is a_11 (no pivoting)
    if a_11 >= α*ω_1
        pivot_size = 1
        push!(pivot, (1,1))
        
    else

        # SECOND SEARCH
        ω_r = max_offdiagonal(A, r)[1]  # returns max element in any off diagonal element in column r


        # The pivot is a_11 (no pivoting)
        if (a_11*ω_r >= α*ω_1^2) 
            pivot_size = 1
            push!(pivot, (1,1))

        # The pivot is a_rr 
        elseif abs(A[r,r]) >= α*ω_r
            pivot_size = 1
            push!(pivot, (1,r))

        # Pivoting between 2 and r (line and column)
        else
            pivot_size = 2
            push!(pivot, (2,r))
        end
    end

    return pivot, pivot_size
end
