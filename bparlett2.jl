using LinearAlgebra
include("max_subdiagonal.jl")

"""
Complete pivoting strategy
"""
function bparlett2(A::Hermitian{T}) where T

    α = (1+sqrt(17))/8
    n = size(A,1)

    pivot = []
    pivot_size = 0

    max_element = abs(A[1,1])
    q = 1
    p = 1

    max_diag = max_element
    r = 1

    # Find max element in sub-diagonal matrix and max diagonal element
    for j in 1:n
        for i in j:n

            # Max on the diagonal : μ_1 = max |a_ij|
            if (i == j)
                if(abs(A[i,j]) > max_diag)
                    max_diag = abs(A[i,j])
                    r = i
                end
            end

            # Max between all elements
            if(abs(A[i,j]) > max_element) 
                max_element = abs(A[i,j])
                p = i
                q = j
            end
        end
    end

    # Pivoting between 1 and r (line and column)
    if (max_diag >= α*max_element)
        pivot_size = 1
        push!(pivot, (1,r))

    # 1) Pivoting between 1 and p (line and column)
    # 2) Pivoting between 2 and q (line and column)
    else
        pivot_size = 2
        push!(pivot, (1,p))
        push!(pivot, (2,q))          
    end

    
    return pivot, pivot_size
end