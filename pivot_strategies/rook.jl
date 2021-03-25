using LinearAlgebra
include("max_subdiagonal.jl")

"""
Rook pivoting strategy
"""
function rook(A::Hermitian{T}) where T
    
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

            ω_r, r = max_subdiagonal(A, i)

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