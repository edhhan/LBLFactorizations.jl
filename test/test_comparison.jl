#Comparing LDLT bkaufmann with Julia bunchkaufmann
using LBLFactorizations
using LinearAlgebra
using Test


#Doesn't work, bunchkaufman doesn't seem to do the same strategy as our
@testset begin

        for _ = 1:5
            for n = 4:100
                             
                A = Hermitian(rand(Float64, n,n).*100, :L)
                F = lbl(A, "bkaufmann")
                FJ = bunchkaufman(A)
                println(norm(F.L - FJ.L))
                println(norm(F.B - FJ.D))

                @test norm(F.L - FJ.L) ≤ sqrt(eps()) * norm(F.L) && norm(F.B - FJ.D) ≤ sqrt(eps()) * norm(F.B)
            end

        end

end



