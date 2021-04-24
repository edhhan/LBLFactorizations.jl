using LinearAlgebra
using ProfileSVG
using CSV
using DataFrames
using Plots
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

#WARMUP
Z = Hermitian(rand(Float64, 2,2).*100)
b=rand(2).*100
F=lbl(Z, "rook")
lbl_solve(F, b)
BK=bunchkaufman(Z)
BK\b

#Test
test_data1 = Float64[]
test_data2 = Float64[]
test_data3 = Float64[]
test_data4 = Float64[]
test_data5 = Float64[]
test_data6 = Float64[]
test_data7 = Float64[]
test_data8 = Float64[]

DIM_START=50
DIM_JUMP=50
DIM_MAX=500

for n=DIM_START:DIM_JUMP:DIM_MAX

    A = Hermitian(rand(Float64, n,n).*100)
    b=rand(n).*100

    t_1 = @elapsed F1=lbl(A, "bkaufmann")
    push!(test_data1, t_1)
    t_2 = @elapsed F2=lbl(A, "bparlett")
    push!(test_data2, t_2)
    t_3 = @elapsed F3=lbl(A, "rook")
    push!(test_data3, t_3)

    t_4 = @elapsed x1=lbl_solve(F1, b)
    push!(test_data4, t_4)
    t_5 = @elapsed x2=lbl_solve(F1, b)
    push!(test_data5, t_5)
    t_6 = @elapsed x3=lbl_solve(F1, b)
    push!(test_data6, t_6)

    #bunchkaufman
    t_7 = @elapsed F4=bunchkaufman(A)
    push!(test_data7, t_7)

    t_8 = @elapsed x4=F4\b
    push!(test_data8, t_8)

    println("Erreur1: ",norm(x1-x4)," ", norm(x2-x4)," ", norm(x3-x4))

    println("lbl bkaufmann : ", t_1)
    println("lbl bparlett: ", t_2)
    println("lbl rook: ", t_3)

    println("lbl bkaufmann solver: ", t_4)
    println("lbl bparlett solver: ", t_5)
    println("lbl rook solver: ", t_6)

    println("Bunchkaufman from Julia construction: ",t_7)
    println("Bunchkaufman from Julia Solver: ",t_8)

    println("Test ", Int64(n/50), " done")


end


    println("Everything done")

    # Creating DataFrame
    data_to_CSV = DataFrame(bkaufmann = test_data1, bparlett = test_data2, rook = test_data3, bkaufmannSolved = test_data4, bparlettSolved = test_data5, rookSolved = test_data6, bunchkaufmanJConstructor = test_data7, bunchkaufmanJSolved= test_data8)  
    filename = "LBL_Tests.csv"

    # modifying the content of "LBL_Tests.csv" using wite method
    CSV.write(filename, data_to_CSV)


    datatest=CSV.File("LBL_Tests.csv")
    df = DataFrame(datatest)

    datatograph1=hcat(df[:,1],df[:,2],df[:,3],df[:,7]*200)
    datatograph2=hcat(df[:,4],df[:,5],df[:,6],df[:,8])
    #Rouler lui
    plot(DIM_START:DIM_JUMP:DIM_MAX,datatograph1, title="Temps de calcul selon la dimension \n (Construction) (bunchkaufman solver*200)", label = ["lbl bkaufmann" "lbl bparlett" "lbl rook" "Bunchkaufman constructor" ] , ylabel="Temps(s)",xlabel="Dimension", legend=:topleft)
    #Rouler lui après
    plot(DIM_START:DIM_JUMP:DIM_MAX,datatograph2, title="Temps de calcul selon la dimension ", label = ["lbl bkaufmann solver" "lbl bparlett solver" "lbl rook solver" "Bunchkaufman solver"] , ylabel="Temps(s)",xlabel="Dimension", legend=:topleft)