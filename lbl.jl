include("pivot_strategies/bkaufmann.jl")
include("pivot_strategies/bparlett.jl")
include("pivot_strategies/rook.jl")
include("LBL_structure.jl")
include("PermuteMatrix.jl")

import Base.push!

using LinearAlgebra
using Random


"""
Wrapper function for pivoting strategy
"""
function pivoting(A::Hermitian{T}, strategy::String) where T

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
"""
function inv_E(E::Union{AbstractMatrix{T}, Array{T}, Float64}, s::Int) where T

    if s==1 || s==0
        return 1/E[1,1]
    elseif s==2
        return inv(E)
    end

end


"""
PAP^T = [E C^* ; C K]

LBL^* Factorization based on 
"""
function lbl(A::Hermitian{T}; strategy::String="rook") where T

    if !ishermitian(A)
        return @error("LBL* factorization only works on hermitian matrices")
    end

    if !(strategy in ["rook", "bparlett", "bkaufmann" ])
        return @error("Invalid pivoting strategy.\nChoose string::strategy ∈ {rook, bparlett, bkaufmann}.")
    end

    # Initialize working matrix
    hat_A = deepcopy(A)

    # Initiliaze data-structure factorization
    n = size(A)[1]
    F = LBL(LowerTriangular{Float64}(zeros(n,n)), zeros(n, n), strategy)

    # Initialize loop variable : undefinite number of iteration
    s = 1

    while s <= n 

        hat_n = size(hat_A)[1]

        # Pivoting
        pivot, pivot_size = pivoting(hat_A, strategy) 

        # Special case for E and L : skip, no permutation matrix required A = [E C^* ; C B]
        if pivot_size == 0 

            # We assign -1 to treat it as a special case
            push_pivot!(F, -1)
            E = hat_A[1,1]
            C = hat_A[2:end, 1]
            B = hat_A[2:end, 2:end]
            L_special_case = vcat(1, zeros(hat_n-1)) #Special case on L

        else

            push_pivot!(F, pivot)

            # If pivot==[(1,1)] then permutation matrix is identity, so we skip that case
            # We direcetly have hat_A = [E C^* ; C B] without any permutations
            if !(pivot == [(1,1)])
                
                # Construct permutation matrix P
                P = Matrix(1.0*I, hat_n, hat_n)

                # pivot is an array of 1 or 2 tuples p
                # p is a tuple of two indices and  
                for p in pivot

                    idx1 = p[1]
                    idx2 = p[2]

                    # Permutations on lines : must be done first [1]
                    temp = P[idx1,:]
                    P[idx1, :] = P[idx2,:]
                    P[idx2, :] = temp

                    # Permutation on columns
                    temp = P[:,idx1]
                    P[:, idx1] = P[:, idx2]
                    P[:, idx2] = temp

                end

                # Apply permuations on working matrix
                hat_A = P*hat_A*P'
               
            end
            
            # With permutation get bloc-matrices from PAP^T = [E C^* ; C B]
            E = hat_A[1:pivot_size, 1:pivot_size]
            C = hat_A[(pivot_size+1):end, 1:pivot_size]
            B = hat_A[(pivot_size+1):end, (pivot_size+1):end]
           
        end

        # Construction of columns of L and B matrices
        E⁻¹ = inv_E(E, pivot_size)
        
        # Special case, where s=1 and no permutation was required
        if pivot_size==0
            F.B[s,s] = E
            F.L[s:end,s] = L_special_case
            #F.L[s:end,s] = vcat(Matrix(1.0*I, 1, 1), C*E⁻¹ )
        else
            # If pivot_size=1, then s+pivot_size-1 = s      =>     s:(s+pivot_size-1) == s:s
            # If pivot_size=2, then s+pivot_size-1 = s+1    =>     s:(s+pivot_size-1) == s:s+1
            F.B[s:(s+pivot_size-1), s:(s+pivot_size-1)] = E
            F.L[s:end, s:(s+pivot_size-1) ] = vcat(Matrix(1.0*I, pivot_size, pivot_size), C*E⁻¹ )
        end
         
        # Schur complement
        if hat_n > 1
            hat_A = Hermitian(B - C*E⁻¹*C')
        end
        
        # Incremental step depends on the size of pivoting
        if pivot_size==1 || pivot_size==0
            s += 1
        elseif pivot_size==2
            s += 2
        end

    end

    return F
end

# Test

#=
using Test
@testset begin

    for _ = 1:4
    
        for n = 4:100
            A = Hermitian(rand(n,n).*100)
            
            F = lbl(A, strategy="rook")
            @test norm(A - F.L*F.B*F.L') ≤ sqrt(eps()) * norm(A)

            #PAPT = PermuteMatrix(A, F.pivot_array)
            #@test norm(PAPT - F.L*F.B*F.L') ≤ sqrt(eps()) * norm(PAPT)
        end
    
    end
end
=#

A_assignement =[1 2 3 4 ; 2 5 13 28 ; 3 13 55 131 ; 4 28 131 270]
A_assignement = Hermitian(A_assignement)
F = lbl(A_assignement, strategy="bkaufmann")
display(F.L)
display(F.B)




#@test norm(A_assignement - F.L*F.B*F.L') ≤ sqrt(eps()) * norm(A)

