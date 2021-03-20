include("bkaufmann.jl")
include("bparlett.jl")
include("rookpivoting.jl")

using LinearAlgebra


"""
Wrapper function 
"""
function pivoting(A::Hermitian{T}, strategy::String) where T

    if strategy == "rook"
        pivot, pivot_size = rook(A)
    elseif strategy == "bparlett"
        pivot, pivot_size = bparlett(A)
    elseif stragy == "bkaufmann"
        pivot, pivot_size = bkaufmann(A)
    end

    return pivot, pivot_size
end


"""
"""
function inv_E(E::AbstractMatrix, s::Int)

    if s==1
        return 1/E
    elseif s==2
        
        # Swap diagonal-element
        temp = E[1,1]  # element A
        E[1,1] = E[2,2]
        E[2,2] = temp

        # Change sign on off diagonal
        E[1,2] = -E[1,2]
        E[2,1] = -E[2,1]
        
        return 1/det(E) * E
    end

end


"""
LBL^* Factorization based on 
"""
function lbl(A::Hermitian{T} ; strategy::String="rook") where T

    if !(strategy in ["rook", "bparlett", "bkaufmann"])
        @error("Invalid pivoting strategy.\nChoose string::strategy ∈ {rook, bparlett, bkaufmann}.")
    end

    # Initialize matrix
    m,n = size(A)
    L = zeros(n,n)
    B = zeros(n,n)

    # Initialize loop variable : undefinite number of iteration
    s = 0

    while s < n

        # Pivoting
        pivot, pivot_size = pivoting(A, strategy) 

        # Special case : skip, no permutation matrix required
        if pivot == 0 
            # TODO
        else

            # Construct permutation matrix
            permutation_matrix = Matrix(1.0*I, n,n)

            # Permutation matrix is identity
            if pivot == [(1,1)]
                continue

            # Permutation on columns and lines are required
            else
                for p in pivot
                    temp = permutation_matrix[:,p]
                end
            end
        end


        # Incremental step depends on the size of pivoting
        if pivot_size==1
            s += 1
        elseif pivot_size==2
            s += 2
        end
    end


    return L,B
end



