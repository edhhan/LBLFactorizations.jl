using LinearAlgebra
include("max_subdiagonal.jl")
include("max_offdiagonal.jl")

"""
Rook pivoting strategy. Provides the pivot according to the rook pivoting strategy (Ashcraft, Grimes, and Lewis). Based on:
    Accuracy and stability of numerical algorithms, Chapter 11, Higham, Nicholas J, 2002, SIAM

    The chapter describe exactly the following elements in the strategy:

    ω_i is the maximimum magnitude of any off-diagonal entry in column r of the previous iteration. Initialized as the maximum magnitude of any subdiagonal entry in column 1.

    r is the row index of the first (subdiagonal) entry of maximum magnitude in column i. 

    ω_r is the maximimum magnitude of any off-diagonal entry in column r.

    α is a parameter derived by minimizing the bound on the element growth.


    This function returns pivot_size=0 and pivot=-1 if there is only nonzero elements under the diagonal element in column 1.
    Otherwise it returns the size of the pivot in pivot_size and the pivot as an array of tuple of the form (k,l) where rows and column k and l have to be swapped. 
"""
#function rook(A::Hermitian{T}) where T
function rook(A::AbstractMatrix{T}) where T

    α = (1+sqrt(17))/8
    n = size(A,1)
    pivot = []
    pivot_size = 0  # pre-initialization required for while-loop

    #Looks for the largest absolute value in the subdiagonal element of column 1
    ω_i, r = max_subdiagonal(A, 1)

    a_11 = abs(A[1,1])
    
    # Early stopping : no pivoting required and trivial solution for the factorization (special case)
    if (ω_i == 0)
        pivot_size = 0
        push!(pivot, -1)
        return pivot, pivot_size
    end

    # The pivot is a_11 (no pivoting)
    if a_11 >= α*ω_i
        pivot_size = 1
        push!(pivot, (1,1))
    else

        i = 1
        # TODO : check special cas -> i == n
        while(pivot_size == 0)

            r = max_subdiagonal(A, i)[2]  # returns index of max subdiagonal element in column i
            ω_r = max_offdiagonal(A, r)[1]  # returns max element in any off diagonal element in column r

            # Pivoting between 1 and r (line and column)
            if abs(A[r,r]) >= α*ω_r 
                pivot_size = 1
                push!(pivot, (1,r))
                
            # 1) Pivoting between 1 and i (line and column)
            # 2) Pivoting between 2 and r (line and column)
            elseif ω_i == ω_r

                pivot_size = 2
                push!(pivot, (1,i))
                push!(pivot, (2,r))
            
            # Rook-displacement at r
            else
                i = r
                ω_i = ω_r
            end

        end

    end

    return pivot, pivot_size

end