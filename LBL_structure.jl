import Base.push!
using LinearAlgebra

abstract type AbstractLBL{T} end

"""
LBL data struture
"""
mutable struct LBL{T} <: AbstractLBL{T} 
    L::UnitLowerTriangular{T}   
    B::Tridiagonal
    strategy::String
    B_inv::Array{Any,1}
    pivot_array::Array{Any,1}
    permutation::Array{Int64,1}
end

"""
Constructor 
"""
function LBL(L::UnitLowerTriangular{T}, B::Tridiagonal{T},
strategy::String="rook", B_inv::Array{Any,1}=Any[], pivot_array::Array{Any,1} = Any[], permutation::Array{Int64,1}= Int64[]) where T   
    return LBL{T}(L, B, strategy, B_inv, pivot_array, permutation)
end

"""
"""
function push_pivot!(A::LBL{T}, pivot::Union{Array{Any,1}, Array{Tuple{Int64, Int64}, 1}, Int64 }) where T
    push!(A.pivot_array, pivot)
end
"""
"""
function push_B_inv!(A::LBL{T}, inv_E::Tuple{Any,Int64}) where T
    push!(A.B_inv, inv_E)
end
"""
"""
function build_matrix(A::LBL{T}) where T
    return A.L * A.B * A.L'
end