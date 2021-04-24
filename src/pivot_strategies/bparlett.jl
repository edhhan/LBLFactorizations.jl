using LinearAlgebra

export bparlett
"""
Complete pivoting strategy in O(n²). Provides the pivot according to the bunch-parlett pivoting strategy.
A search is done in the subdiagonal part of a matrix A. The pivots on the main diagonal are prioritized.
    
    The chapter 11 in [1] describe the following elements in the strategy:
        -μ_0 : max element in upper triangular with row index p and column index q
        -μ_1 : max diagonal element  with row and column index r
        -α   :  a parameter derived by minimizing the bound on the element growth [1]
    
    Input :
        -A : a matrix 
    Output :
        -pivot_size: size of the pivot, i.e size ∈ {1, 2}
        -pivot : the indices linked to the permutation (pivoting) 
                -if size==1 then pivot = tuple(idx1,idx2)
                -if size==2 then pivot = [tuple(idx1,idx2), tuple(idx3,idx4)] 


    This function returns the size of the pivot in pivot_size and the pivot as an array of tuple of the form (k,l) where rows and column k and l have to be swapped. 
        
Based on :     
[1] N. J. Higham, Accuracy and stability of numerical algorithms (Chapter 11). SIAM, 2002.
"""
function bparlett(A::AbstractMatrix{T}) where T
    α = (1+sqrt(17))/8
    n = size(A,1)

    pivot = []
    pivot_size = 0

    μ_0 = 0
    μ_1 = 0
    r = 1
    p = 1
    q = 1

    # Find max element in sub-diagonal matrix and max diagonal element
    for j in 1:n # columns
        for i in 1:j # lines

            # Max diagonal element 
            if (i == j)
                if(abs(A[i,j]) > μ_1)
                    μ_1 = abs(A[i,j])
                    r = i
                end
            end

            # Max element
            if (abs(A[i,j]) > μ_0) 
                μ_0 = abs(A[i,j])
                p = i
                q = j
            end

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

    
    return pivot, pivot_size
end