using LinearAlgebra
include("bparlett.jl")
include("bkaufmann.jl")

n=5

A=rand(n,n)*1000+rand(n,n)+rand(n,n)+rand(n,n)+rand(n,n)+rand(n,n)-500*I
B=Symmetric(A)

pivots, dimension = bparlett(B)
pivots, dimension = bkaufmann(B)

display(B)
print(pivots,dimension)