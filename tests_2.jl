using LinearAlgebra
using Test

include("pivot_strategies/bparlett.jl")
include("pivot_strategies/bkaufmann.jl")
include("pivot_strategies/rook.jl")
include("lbl.jl")


# Test

#=
@testset begin

    for _ = 1:120
    
        for n = 5:5
            A = Hermitian(rand(n,n).*100)
            
            F = lbl(A, strategy="rook")
            
            #@test norm(A - F.L*F.B*F.L') ≤ sqrt(eps()) * norm(A)

            PAPT = permute_matrix_lbl(A, F.pivot_array)
            @test norm(PAPT - F.L*F.B*F.L') ≤ sqrt(eps()) * norm(PAPT)
        end
    
    end
end
=#


#=
A_assignement =[1 2 3 4 ; 2 5 13 28 ; 3 13 55 131 ; 4 28 131 270]
A_assignement = Hermitian(A_assignement)
F = lbl(A_assignement, strategy="bparlett")
display(F.L)
display(F.B)
@test norm(A_assignement - F.L*F.B*F.L') ≤ sqrt(eps()) * norm(A_assignement)

=#
# Bunch-parlet 
#=function permutation_matrix(pivot,n)

    P=Matrix(1.0I,n,n)
    P[pivot[1], pivot[1]]=0
    P[pivot[2], pivot[2]]=0

    P[pivot[2], pivot[1]]=1
    P[pivot[1], pivot[2]]=1
    return P

end
=#
# Itération 1 
n=100
A = rand(n,n)*1000+rand(n,n)+rand(n,n)+rand(n,n)+rand(n,n)+rand(n,n)-500*I
#A =[1 2 3 4 ; 2 5 13 28 ; 3 13 55 131 ; 4 28 131 270]
display(A)
A = Hermitian(A,:L)

F1 = lbl(A, "bparlett")
F = bunchkaufman(A);


for permutation in F1.permutation_array
    P = permutation_matrix(permutation, n)
    global A = P*A*P'
end

println( norm(A-  build_matrix(F1))  )