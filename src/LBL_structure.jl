import Base.push!
using LinearAlgebra

abstract type AbstractLBL{T} end

"""
LBL data struture that contains the information of a LBL^* factorization (see lbl.jl) and utility information.

Main attributes:
    -L : the triangular matrix of the factorzation (also implicitly representes the L^* matrix of the factorization)
    -B : the bloc-diagonal matrix, i.e. each bloc is 1x1 and 2x2, implemented as a Tridiagonal matrix
Utility attributes:
    -strategy : a pivoting strategy chosen by the user, i.e strategy ∈ {rook, bparlett, bkaufmann}
    -B_inv : the inverse of the block-matrices in B that is used in solving solve(LBL^,b)
    -pivot_array : a array containing the pivot done in the factorization, containing three type of elements : 
                    1) a tuple represents a single pivot for the given iteration (index in pivot_array)
                    2) [tuple, tuple] represents two pivots for the given iteration
                    3) -1 special case where the column L contains 0 in its subdiagonal (bkaufmann and rook)   
    -permutation : indices that describes the permutations done in the factorization : can be used to reconstruct A

"""
mutable struct LBL{T} <: AbstractLBL{T} 
    L::UnitLowerTriangular{T}   
    B::Tridiagonal{T}
    strategy::String
    B_inv::Array{Any,1}
    pivot_array::Array{Any,1}
    permutation::Array{Int64,1}
end


"""
Constructor by parameters that initializes the empty arrays for B_inv, pivot_array and permutation
"""
function LBL(L::UnitLowerTriangular{T}, B::Tridiagonal{T}, strategy::String="rook", B_inv::Array{Any,1}=Any[],
                                     pivot_array::Array{Any,1} = Any[], permutation::Array{Int64,1}=Int64[]) where T   
    return LBL{T}(L, B, strategy, B_inv, pivot_array, permutation)
end

"""
Overload of the push function for the array-attribute pivot_array
"""
function push_pivot!(A::LBL{T}, pivot::Union{Array{Any,1}, Array{Tuple{Int64, Int64}, 1}, Int64}) where T
    push!(A.pivot_array, pivot)
end

"""
Overload of the push function for the array-attribute push_B_inv
"""
function push_B_inv!(A::LBL{T}, inv_E::Any) where T
    push!(A.B_inv, inv_E)
end

