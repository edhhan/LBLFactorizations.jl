include("pivot_strategies/bkaufmann.jl")
include("pivot_strategies/bparlett.jl")
include("pivot_strategies/rook.jl")
include("LBL_structure.jl")

using LinearAlgebra

"""
Return the permutation matrice associated with the pivot
"""
function permutation_matrix(permutation,n)
    P = zeros(n,n)
    for i=1:n
        P[i,permutation[i]]=1
    end

    return P
end
"""
Wrapper function for pivoting strategy
"""
#function pivoting(A::Hermitian{T}, strategy::String) where T
function pivoting(A::AbstractMatrix{T}, strategy::String) where T
    if strategy == "rook"
        pivot, pivot_size = rook(A)
    elseif strategy == "bparlett"
        pivot, pivot_size = bparlett(A)
    elseif strategy == "bkaufmann"
        pivot, pivot_size = bkaufmann(A)
    end

    return pivot, pivot_size
end



"""
LBL^* Factorization based on 
"""
function lbl(A::Hermitian{T}, strategy::String="bkaufmann") where T   
    if !ishermitian(A)
        return @error("LBL* factorization only works on hermitian matrices")
    end

    if !(strategy in ["rook", "bparlett", "bkaufmann" ])
        return @error("Invalid pivoting strategy.\nChoose string::strategy ∈ {rook, bparlett, bkaufmann}.")
    end


    # Initialization
    n = size(A,1)
    hat_n = n
    F = LBL(UnitLowerTriangular(zeros(n,n)), Tridiagonal(zeros(n, n)), strategy)
    #F = LBL(zeros(n,n), zeros(n, n), strategy)
    F.permutation=1:n
    A_prime = Matrix(A) # A_prime cannot be hermitian because of the inplace permutations below

    s = Int64(1)
    while(s < n)

        pivot, pivot_size = pivoting(A_prime, strategy)
        push_pivot!(F, pivot)

        if pivot_size != 0

            for p in pivot 

                if p != (1,1)
                    
                    # Permutation inplace on lines 
                    temp = A_prime[p[1], :]
                    A_prime[p[1], :] = A_prime[p[2], :]
                    A_prime[p[2], :] = temp
                    
                    # Permutation inplace on columns
                    temp = A_prime[:, p[1]]
                    A_prime[:, p[1]] = A_prime[:,p[2]]
                    A_prime[:, p[2]] = temp

                    # Permutation of UnitLowerTriangular matrix L
                    temp = F.L[s+p[1]-1, 1:s+p[1]-2]
                    F.L[s+p[1]-1, 1:s+p[1]-2] = F.L[s+p[2]-1, 1:1:s+p[1]-2] 
                    F.L[s+p[2]-1, 1:s+p[1]-2] = temp 

                    # Updating permutation vector
                    temp=F.permutation[s+p[1]-1]
                    F.permutation[s+p[1]-1]=F.permutation[s+p[2]-1]
                    F.permutation[s+p[2]-1]=temp

                end
            end

        end

        # PAP^T = [E C^* ; C K]
        B = A_prime[(pivot_size+1):hat_n,(pivot_size+1):hat_n] #B
        C = A_prime[(pivot_size+1):hat_n,1:pivot_size]         #C
        E = A_prime[1:pivot_size,1:pivot_size]                 #E

        # Schur complement 
        #A_prime = Hermitian(B - C*inv(E)*C')
        inv_E=inv(A_prime[1:pivot_size,1:pivot_size])
        A_prime = (B - C*inv_E*C') #Prend bcp de temps
        hat_n = size(A_prime,1)

        # Fill factorization columns in L and block-diagonal (1x1 or 2x2) in E
        F.L[s+pivot_size:end,s:s+pivot_size-1] = C*inv_E
        F.B[s:s+pivot_size-1,s:s+pivot_size-1] = E
        push_B_inv!(F, (inv_E,s))
      
        # Incremental step depends on the size of pivoting
        if pivot_size==1 || pivot_size==0
            s += 1
        elseif pivot_size==2
            s += 2
        end

    end
    
    # Last step
    if s == n
        F.L[n, n] = 1
        F.B[n, n] = A_prime[1,1]
        inv_E=inv(A_prime[1,1])
        push_B_inv!(F, (inv_E,s))
    end

    return F
end
