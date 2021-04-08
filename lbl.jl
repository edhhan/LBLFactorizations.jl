include("pivot_strategies/bkaufmann.jl")
include("pivot_strategies/bparlett.jl")
include("pivot_strategies/rook.jl")
include("LBL_structure.jl")

using LinearAlgebra


function permutation_matrix(pivot,n)
    P = Matrix(1.0I,n,n)
    P[pivot[1], pivot[1]] = 0.0
    P[pivot[2], pivot[2]] = 0.0
    P[pivot[2], pivot[1]] = 1.0
    P[pivot[1], pivot[2]] = 1.0
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
#function lbl(A::Hermitian{T}, strategy::String="bkaufmann") where T
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
    #F = LBL(UnitLowerTriangular(zeros(n,n)), zeros(n, n), strategy)
    F = LBL(zeros(n,n), zeros(n, n), strategy)
    A_prime = Matrix(A) # A_prime cannot be hermitian because of the inplace permutations below

    s = 1
    while(s < n)

        pivot, pivot_size = pivoting(A_prime, strategy)
        push_pivot!(F, pivot)

        if pivot_size != 0

            for p in pivot

                if p != (1,1)

                    #A_prime = P*A_prime*P'
                    
                    # Permutation inplace on lines
                    temp = A_prime[p[1], :]
                    A_prime[p[1], :] = A_prime[p[2], :]
                    A_prime[p[2], :] = temp
                    
                    # Permutation inplace on columns
                    temp = A_prime[:, p[1]]
                    A_prime[:, p[1]] = A_prime[:,p[2]]
                    A_prime[:, p[2]] = temp

                    #P_augmented = Matrix(1.0*I, n,n)
                    #P_augmented[s:end,s:end] = P
                    temp = F.L[s+p[1]-1, :]
                    F.L[s+p[1]-1, :] = F.L[s+p[2]-1, :]
                    F.L[s+p[2]-1, :] = temp

                    
                    #F.L = UnitLowerTriangular(P_augmented*F.L) # P*hat_A : permute lines
                    push_permutation!(F, (s+p[1]-1, s+p[2]-1) ) 
                    

                end
            end

        end

        # PAP^T = [E C^* ; C K]
        B = A_prime[(pivot_size+1):end,(pivot_size+1):end]
        C = A_prime[(pivot_size+1):end,1:pivot_size]
        E = A_prime[1:pivot_size,1:pivot_size]

        # Schur complement 
        #A_prime = Hermitian(B - C*inv(E)*C')
        A_prime = (B - C*inv(E)*C')
        hat_n = size(A_prime,1)

        # Fill factorization columns
        F.L[s:end,s:s+pivot_size-1] = vcat(Matrix(I, pivot_size, pivot_size), C*inv(E) )
        F.B[s:s+pivot_size-1,s:s+pivot_size-1] = E
        
      
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
    end

    return F
end


export ldiv!
function lbl_solve!(LBL::AbstractLBL, b:;AbstractVector) where T   
    

end