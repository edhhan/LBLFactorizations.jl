using LinearAlgebra
using Test
using BenchmarkTools
using ProfileSVG

include("../lbl.jl")
include("../lbl_solve.jl")

n=10000
A = Hermitian(rand(Float64, n,n).*100)
ProfileSVG.@profview  F = lbl(A, "rook")
b=rand(n).*100


#ProfileSVG.@profview lbl_solve(F, b)
println("DONE")
#ProfileSVG.@profview lbl_solve(F, b)
#time1 = @elapsed x1=lbl_solve(F, b)
#time2= @elapsed x2=A\b

#display(x1)
#display(x2)
#println("Erreur",norm(x1-x2))
#println("Temps lblsolve: ",time1)
#println("Temps julia: ",time2)
