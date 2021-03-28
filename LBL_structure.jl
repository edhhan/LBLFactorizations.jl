import Base.push!
using LinearAlgebra


abstract type AbstractLBL{T} end

"""
LBL data struture
"""
mutable struct LBL{T} <: AbstractLBL{T} 
    #L::LowerTriangular{T}
    L::AbstractMatrix{T}
    B::AbstractMatrix{T}
    
    strategy::String
    pivot_array::Array{Any,1}
    permutation_array::Array{Any,1}
end

"""
Constructor 
"""
#function LBL(L::LowerTriangular{T}, B::AbstractMatrix{T}, strategy::String="rook", pivot_array::Array{Any,1} = Any[] ) where T 
function LBL(L::AbstractMatrix{T}, B::AbstractMatrix{T},
             strategy::String="rook", pivot_array::Array{Any,1} = Any[], permutation_array::Array{Any,1} = Any[]  ) where T   
    return LBL{T}(L, B, strategy, pivot_array, permutation_array )
end


"""
"""
function push_pivot!(A::LBL{T}, pivot::Union{Array{Any,1}, Array{Tuple{Int64, Int64}, 1}, Int64 }) where T
    push!(A.pivot_array, pivot)
end

"""
"""
function push_permutation!(A::LBL{T}, pivot::Tuple{Int64, Int64}) where T
    push!(A.permutation_array, pivot)
end

"""
"""
function build_matrix(A::LBL{T}) where T
    return A.L * A.B * A.L'
end