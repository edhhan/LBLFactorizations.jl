module LBLFactorizations

export lbl, lbl_solve, bkaufmann, bparlett, rook, build_matrix, permutation_matrix

include("pivot_strategies/bkaufmann.jl")
include("pivot_strategies/bparlett.jl")
include("pivot_strategies/rook.jl")
include("LBL_structure.jl")
include("lbl_solve.jl")

using LinearAlgebra

"""
A LBL^* factorization, also known as a LDL^* bloc-factorization, where the L matrix is UnitLowerTriangular
and the B matrix is a bloc-diagonal matrix, i.e. each block is 1x1 or 2x2.

The LBL^* factorization is a generalization of the common LDL^T factorization : it handles hermitian 
and also indefinite matrices.

The factorization requires pivoting. Three strategies are proposed:
    -bparlett : a pivoting strategy with thorough search that is expensive : search in the subdiagonal O(n²)
                 (see pivot_strategies/bparlett.jl for more details)
    
    -bkaufmann : a partial pivoting strategy that is rapid and less precise : search only two columns O(n)
                 (see pivot_strategies/bkaufmann.jl for more details)

    -rook : a compromise between the two strategies, that has a complexity between O(n) and O(n²) [2]

Each strategy is relevant : it depends on the the structure of the matrix and the user needs.


    Input :
        -A : a hermitian matrix
        -strategy : a pivoting strategy chosen by the user, i.e. strategy ∈ {rook, bparlett, bkaufmann}. 
                    Default strategy is set as rook (compromise).

    Output :
        -F : a LBL data structure containing multiple attributes (see LBL_structure.jl for more details)


[1] N. J. Higham, Accuracy and stability of numerical algorithms (Chapter 11). SIAM, 2002.

[2] G. Poole and L. Neal, “The rook’s pivoting strategy,”Journal of Computational and Applied Mathematics,
    vol. 123, no. 1-2, pp. 353–369, 2000.
"""
function lbl(A::Union{Hermitian{T}, AbstractMatrix{T}}, strategy::String="rook") where T   

    if !ishermitian(A)
        return @error("LBL* factorization only works on hermitian matrices")
    end

    if !(strategy in ["rook", "bparlett", "bkaufmann" ])
        return @error("Invalid pivoting strategy.\nChoose string::strategy ∈ {rook, bparlett, bkaufmann}.")
    end


    # Initialization
    n = size(A,1)
    hat_n = n
    Zeros=zeros(n, n)
    F = LBL(UnitLowerTriangular(Zeros), Tridiagonal(Zeros), strategy)
    F.permutation = 1:n

    A_prime = Matrix(A) # A_prime cannot be hermitian because of the inplace permutations below

    s = 1
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

        # PAP^T = [E C^* 
        #          C B ]
        # Matrices below are identified in the PAP^T matrix 
        B = A_prime[(pivot_size+1):hat_n, (pivot_size+1):hat_n] 
        C = A_prime[(pivot_size+1):hat_n, 1:pivot_size]         
        E = A_prime[1:pivot_size, 1:pivot_size]                 

        # Schur complement 
        inv_E = inv(A_prime[1:pivot_size,1:pivot_size])
        A_prime = B - C*inv_E*C' # NOTE : botte-neck in the FlameGraph
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

"""
Utility function that constructs the permutation matrix associated to the a vector of permutation.
The function is used in the tests/test_assignement.jl file in which we reconstruct the matrix A, i.e A=L*B*L'

    Input :
        -permutation : a vector of permutation (permutation attribute of the LBL_structure), i.e. LBL.permutation
        -n : size of the square matrix

    Output :
        -P : a permutation matrix
"""
function permutation_matrix(permutation, n)
    
    P = zeros(n,n)
    for i = 1:n
        P[i,permutation[i]] = 1
    end

    return P
end


"""
Wrapper function pivoting strategy in lbl() : enhances the readability

    Input : 
        -A : an AbstractMatrix, called A_prime in the lbl() function
        -strategy : a pivoting strategy, i.e. strategy ∈ {rook, bparlett, bkaufmann}
                    (see pivot_strategies dir for more info)

    Output :
        -pivot : a tuple (idx1, idx2) that describes the indices to permute
        -pivot_size : a size, i.e size ∈ {0, 1, 2}, that enables us to treat different cases in the lbl() function
"""
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
Utility function for tests 
"""
function build_matrix(A::LBL{T}) where T
    return A.L * A.B * A.L'
end


# Second version of lbl to test modifications : robust Hermitian type for A_prime (less effective/slower)
function lbl_v2(A::Union{Hermitian{T}, AbstractMatrix{T}}, strategy::String="rook") where T   

    if !ishermitian(A)
        return @error("LBL* factorization only works on hermitian matrices")
    end

    if !(strategy in ["rook", "bparlett", "bkaufmann" ])
        return @error("Invalid pivoting strategy.\nChoose string::strategy ∈ {rook, bparlett, bkaufmann}.")
    end


    # Initialization
    n = size(A,1)
    hat_n = n
    Zeros=zeros(n, n)
    F = LBL(UnitLowerTriangular(Zeros), Tridiagonal(Zeros), strategy)
    F.permutation = 1:n


    if A.uplo=='L'
        A_prime=deepcopy(Hermitian(A.data',:U))
    else

    A_prime = deepcopy(A)
    end
    s = 1
    while(s < n)

        pivot, pivot_size = pivoting(A_prime, strategy)
        push_pivot!(F, pivot)
        
        if pivot_size != 0
            for p in pivot 
                if p != (1,1)

                    p1=p[1]
                    p2=p[2]

                    temp=A_prime.data[p1,p1+1:p2-1]
                    A_prime.data[p1,p1+1:p2-1]=A_prime.data[p1+1:p2-1,p2]'
                    A_prime.data[p1+1:p2-1,p2]=temp'
                            
                    A_prime.data[p1,p2]=A_prime.data[p1,p2]'
                            
                    temp=A_prime.data[p1,p2+1:end]
                    A_prime.data[p1,p2+1:end]=A_prime.data[p2,p2+1:end]
                    A_prime.data[p2,p2+1:end]=temp
                        
                    temp=A_prime.data[1:p1-1,p1]
                    A_prime.data[1:p1-1,p1]=A_prime.data[1:p1-1,p2]
                    A_prime.data[1:p1-1,p2]=temp   
                        
                    temp=A_prime.data[p1,p1]
                    A_prime.data[p1,p1]=A_prime.data[p2,p2]
                    A_prime.data[p2,p2]=temp

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

        # PAP^T = [E C^* 
        #          C B ]
        # Matrices below are identified in the PAP^T matrix 
        B = Hermitian(A_prime[(pivot_size+1):hat_n, (pivot_size+1):hat_n]) 
        C = A_prime[(pivot_size+1):hat_n, 1:pivot_size]   
        E = Hermitian(A_prime[1:pivot_size, 1:pivot_size])                 

        # Schur complement 
        inv_E = inv(E)
        A_prime = Hermitian(B - C*inv_E*C', :U) # NOTE : botte-neck in the FlameGraph
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

end # module
