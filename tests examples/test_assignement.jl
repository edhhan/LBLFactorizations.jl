using LinearAlgebra
using Test
using LBLFactorizations

@testset begin

    for strategy in ["bparlett", "bkaufmann", "rook"]

        A = Hermitian([1 2 3 4 ; 2 5 13 28 ; 3 13 55 131 ; 4 28 131 270])
        n = size(A)[1]

        F = lbl(A, strategy)
        F_build = build_matrix(F)

        # Permute matrix A with permutations within LBL^* factorization
        P = permutation_matrix(F.permutation, n)
        A = P*A*P'
        
        @test norm(A - F_build) ≤ sqrt(eps()) * norm(A)
        display(F_build)
    end

end