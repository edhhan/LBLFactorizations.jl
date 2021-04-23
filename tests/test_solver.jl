using LinearAlgebra
using ProfileSVG
#using ProfileView

include("../lbl.jl")
include("../lbl_solve.jl")


n=1000 #Choisir dimension du probleme, MAX 2000, trop de memoire requise WORK IN PROGRESS

     
#=
    A = Hermitian(rand(Float64, n,n).*100)
    Z = Hermitian(rand(Float64, n,n).*100)

    @profview lbl(A, "rook")
    #ProfileSVG.@profview lbl(A, "rook")
=#
#=
    A = Hermitian(rand(Float64, n,n).*100)
    F=lbl(A, "rook")
    b=rand(n).*100
    ProfileSVG.@profview lbl_solve(F, b)
=#


    A = Hermitian(rand(Float64, n,n).*100)
    Z = Hermitian(rand(Float64, 2,2).*100)
    b=rand(n).*100

    time0=@elapsed F2=lbl(Z, "rook")
    time0 += @elapsed F1=lbl(A, "rook")

    time2= @elapsed lbl_solve(F2, b)
    time2 += @elapsed x1=lbl_solve(F1, b)

    time4 = @elapsed x3=A\b

    println("Erreur1: ",norm(x1-x3))

    println("Temps lbl: ",time0)

    println("Temps lbl_solve: ",time2)

    println("Temps julia: ",time4)

    println("DONE")

