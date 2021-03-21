import Base.push!
using LinearAlgebra


abstract type AbstractLBL{T} end

"""
LBL data struture
"""
mutable struct LBL{T} <: AbstractLBL{T} 
    L::LowerTriangular{T}
    B::AbstractMatrix{T}
    strategy::String
    pivot_array::Array{Any,1}
end

"""
Constructor 
"""
function LBL(L::LowerTriangular{T}, B::AbstractMatrix{T}, strategy::String="rook", pivot_array::Array{Any,1} = Any[] ) where T 
    return LBL{T}(L, B, strategy, pivot_array )
end


"""
"""
function push_pivot!(A::LBL{T}, pivot::Union{Array{Any,1}, Array{Tuple{Int64, Int64}, 1}, Int64 }) where T
    push!(A.pivot_array, pivot)
end