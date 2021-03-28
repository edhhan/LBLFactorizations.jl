#Comparing LDLT bkaufmann with Julia bunchkaufmann
include("../pivot_strategies/bparlett.jl")
include("../pivot_strategies/bkaufmann.jl")
include("../pivot_strategies/rook.jl")
include("../lbl.jl")
using LinearAlgebra

#@testset begin

        for _ = 1:5
            for n = 4:9
                             
                A = Hermitian(rand(Float64, n,n).*100, :L)
                F = lbl(A, "bkaufmann")
                FJ = bunchkaufman(A)
                println(norm(F.L - FJ.L))
                println(norm(F.B - FJ.D))

                #@test norm(F.L - FJ.L) ≤ sqrt(eps()) * norm(F.L) && norm(F.B - FJ.D) ≤ sqrt(eps()) * norm(F.B)
            end

        end

#end



