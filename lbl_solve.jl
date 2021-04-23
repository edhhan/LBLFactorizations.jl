include("LBL_structure.jl")

using LinearAlgebra
#Solving the Bx=b system using precalculated inverses of the blocks of B, B_inv[i][2] contains the index of
# position of the block inverses in the block matrix B.d
#To parallelize
function solve_block_diagonal(B_inv::Array{Any,1}, b::AbstractVector)
    n=length(B_inv)
    for i=1:n
        dim=size(B_inv[i][1],1)
        #Multiplying the member for the right (b) by the inverse of the blocks of B
        b[B_inv[i][2]:B_inv[i][2]+dim-1]=B_inv[i][1]*b[B_inv[i][2]:B_inv[i][2]+dim-1]
    end
    return b
end


function lbl_solve(LBL::AbstractLBL, b::AbstractVector)

    b=b[LBL.permutation] 
    #Triangular Descent
    z = LBL.L \ b 

    #Independant bloc diagonal
    y = solve_block_diagonal(LBL.B_inv,z) #LBL.B \ z 

    #Triangular Hike
    x = LBL.L' \ y 
    x = invpermute!(x,LBL.permutation)
    return x
end