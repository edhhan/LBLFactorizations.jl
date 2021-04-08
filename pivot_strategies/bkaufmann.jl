using LinearAlgebra
include("max_subdiagonal.jl")

"""
Partial pivoting strategy. Provides the pivot according to the bunch-kaufmann pivoting strategy. Based on:
Accuracy and stability of numerical algorithms, Chapter 11, Higham, Nicholas J, 2002, SIAM

The chapter describe the following elements in the strategy:

ω_1 is the maximimum magnitude of any subdiagonal entry in column 1

r is the row index of ω_1
ω_r is the maximimum magnitude of any off-diagonal entry in column r

α is a parameter derived by minimizing the bound on the element growth.

This function returns pivot_size=0 and pivot=-1 if there is only nonzero elements under the diagonal element in column 1.
Otherwise it returns the size of the pivot in pivot_size and the pivot as an array of tuple of the form (k,l) where rows and column k and l have to be swapped. 
"""
#function bkaufmann(A::Hermitian{T}) where T
function bkaufmann(A::AbstractMatrix{T}) where T
    α = (1+sqrt(17))/8
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
