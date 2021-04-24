using LinearAlgebra
using ProfileSVG
using ProfileView

include("../lbl.jl")
include("../lbl_solve.jl")

#Choose the dimension of the problem to run get a profile from. As of now, n shouldn't be greater than 2000 because lbl takes too much time to factorize it (may freeze computer)
n=1000 #MAX 2000

A = Hermitian(rand(Float64, n,n).*100)

ProfileView.@profview lbl(A, "rook")
#ProfileSVG.@profview lbl(A, "rook")

#=
    A = Hermitian(rand(Float64, n,n).*100)
    F=lbl(A, "rook")
    b=rand(n).*100
    ProfileSVG.@profview lbl_solve(F, b)
=#