using LinearAlgebra
include("max_subdiagonal.jl")

"""
Complete pivoting strategy
"""
function bparlett(A::Hermitian{T}) where T

    α=(1+sqrt(17))/8
    n=size(A,1)

    μ_0 = 0
    μ_1 = 0

    pivot = []
    pivot_size = 0

    for j in 1:n
        for i in j:n

            v = abs(A[i,j])

            # Max on the diagonal : μ_1 = max |a_ij|
            if (i == j)
                if(v > μ_1)
                    μ_1 = A[i,j]
                    r = i
                end
            end
            # Max between all elements
            if(v > μ_0) 
                μ_0 = A[i,j]
            end
        end
        

        # Pivoting between 1 and r (line and column)
        if (μ_1 >= α*μ_0)
            pivot_size = 1
            push!(pivot, (1,r))

        # 1) Pivoting between 1 and p (line and column)
        # 2) Pivoting between 2 and q (line and column)
        else
            pivot_size = 2
            push!(pivot, (1,p))
            push!(pivot, (2,q))          
        end

    end

    return pivot, pivot_size
end