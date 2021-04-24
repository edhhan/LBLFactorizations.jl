using LinearAlgebra
using Test

include("../pivot_strategies/bparlett.jl")
include("../pivot_strategies/bkaufmann.jl")
include("../pivot_strategies/rook.jl")
include("../lbl.jl")

@testset begin 

    for strategy in ["bparlett", "bkaufmann", "rook"]

        for _ = 1:5
            for n = 4:100 
                             
                A = Hermitian(rand(Float64, n,n).*100, :L)
                F = lbl(A, strategy)
                F_build = build_matrix(F)

                # Permute matrix A with permutations within LBL^* factorization
                P = permutation_matrix(F.permutation, n)
                A = P*A*P'
                
                @test norm(A - F_build) ≤ sqrt(eps()) * norm(A)
            end

        end

    end

end
