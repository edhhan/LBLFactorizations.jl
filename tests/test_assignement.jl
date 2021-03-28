using LinearAlgebra
using Test

include("../pivot_strategies/bparlett.jl")
include("../pivot_strategies/bkaufmann.jl")
include("../pivot_strategies/rook.jl")
include("../lbl.jl")

@testset begin

    for strategy in ["bparlett", "bkaufmann", "rook"]

        A = Hermitian([1 2 3 4 ; 2 5 13 28 ; 3 13 55 131 ; 4 28 131 270])
        n = size(A)[1]

        F = lbl(A, strategy)
        F_build = build_matrix(F)

        # Permute matrix A with permutations within LBL^* factorization
        for permutation in F.permutation_array
            P = permutation_matrix(permutation, n)
            A = P*A*P'
        end
        
        @test norm(A - F_build) ≤ sqrt(eps()) * norm(A)
        display(F_build)
    end

end