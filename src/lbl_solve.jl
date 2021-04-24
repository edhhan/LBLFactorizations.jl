include("LBL_structure.jl")

using LinearAlgebra

"""
Solver for a (LBL^*)x = b system using a LBL data structure where a LBL^* factorization has been applied (see lbl.jl).


    Input:
        -LBL : a LBL structure containing a factorized matrix 
        -b : RHS vector 

    Output:
        -x : solved

    Approach : the (LBL^*) x = b  is divided into three subproblems :
        1) Solve for z : a triangular Descent with Lz=B
        2) Solve for y : a independant bloc-diagonal By=Z
        3) Solve for x : a triangular hike L^* x = y

    Thus, LBL^*x = LB(y) = L(z) = b

"""
function lbl_solve(LBL::AbstractLBL, b::AbstractVector)

    b = b[LBL.permutation] 

    # 1) Triangular Descent
    z = LBL.L \ b 

    # 2) Independant bloc diagonal : possibility to be parallelized
    y = solve_block_diagonal!(LBL.B_inv,z) #LBL.B \ z 

    # 3) Triangular Hike
    x = LBL.L' \ y 
    x = invpermute!(x,LBL.permutation)

    return x
end



"""
Sub-solver for a Bx=b system using precalculated inverses of the blocks of B, where B is block-diagonal
with each matrix 2x2 or 1x1. B_inv[i][2] contains the index of position of the block inverses in the block matrix B.d

    Input:
        -B_inv : array the inverses of the block-matrix (stored in LBL.B_inv) 
        -b : RHS vector solved inplace

    Output:
        -b : b is solved inplace, thus, as an output, it represents x in Bx=b 

TO DO : parallelize, @avxt thread=x? with LoopVectorization.jl
"""
function solve_block_diagonal!(B_inv::Array{Any,1}, b::AbstractVector)
    n = length(B_inv)
    for i = 1:n
        dim = size(B_inv[i][1],1)
        
        #Multiplying the member for the right (b) by the inverse of the blocks of B
        b[B_inv[i][2]:B_inv[i][2]+dim-1] = B_inv[i][1]*b[B_inv[i][2]:B_inv[i][2]+dim-1]
    end
    return b
end
