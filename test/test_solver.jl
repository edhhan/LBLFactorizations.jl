using LinearAlgebra
using ProfileSVG
using CSV
using DataFrames
using Plots
using ProfileView


#=The testing loop below gets the time required to the lbl factorization for:
- Factorizing A with the bunch-kaufmann strategy from lbl
- Factorizing A with the bunch-parlett strategy from lbl
- Factorizing A with the rook strategy from lbl
- Factorizing A with the bunchkaufman function from LinearAlgebra

- Solving Ax=b with the factored A via the bunch-kaufmann strategy from lbl and with solver from lbl
- Solving Ax=b with the factored A via the bunch-parlett strategy from lbl and with solver from lbl
- Solving Ax=b with the factored A via the rook strategy from lbl and with solver from lbl
- Solving Ax=b with the factored A via the bunchkaufman function from LinearAlgebra and its solver
=#

#=
#WARMUP
Z = Hermitian(rand(Float64, 2,2).*100)
b=rand(2).*100
F=lbl(Z, "rook")
lbl_solve(F, b)
BK=bunchkaufman(Z)
BK\b

#Test data, to create CSV down below
test_data9 = Int64[]
test_data1 = Float64[]
test_data2 = Float64[]
test_data3 = Float64[]
test_data4 = Float64[]
test_data5 = Float64[]
test_data6 = Float64[]
test_data7 = Float64[]
test_data8 = Float64[]

#Sampling options
#DIM_JUMP is the period of the Sampling
#DIM_MAX is the maximum possible dimension of the considered problems. As of now, it shouldn't be greater than 2000 because lbl takes too much time to factorize it (may freeze computer)
#DIM_START is the starting dimension
DIM_START=50
DIM_JUMP=50
DIM_MAX=100
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
    push!(test_data9, n)

    println("Erreur1: ",norm(x1-x4)," ", norm(x2-x4)," ", norm(x3-x4))

    println("lbl bkaufmann : ", t_1)
    println("lbl bparlett: ", t_2)
    println("lbl rook: ", t_3)

    println("lbl bkaufmann solver: ", t_4)
    println("lbl bparlett solver: ", t_5)
    println("lbl rook solver: ", t_6)

    println("Bunchkaufman from Julia construction: ",t_7)
    println("Bunchkaufman from Julia Solver: ",t_8)

    println("Test ", n, " done")


end
println("Everything done")
=#
#Uncomment to export to CSV and plot graphs


#=
    # Creating DataFrame
    data_to_CSV = DataFrame(bkaufmann = test_data1, bparlett = test_data2, rook = test_data3, bkaufmannSolved = test_data4, bparlettSolved = test_data5, rookSolved = test_data6, bunchkaufmanJConstructor = test_data7, bunchkaufmanJSolved= test_data8, dim=test_data9) 
    # Name of the file to export to 
    filename = "LBL_Tests.csv"

    # Writing the data in CSV format
    CSV.write(filename, data_to_CSV)
=#
    # Collecting the data in CSV format, should be in the same directory
    datatest=CSV.File("LBL_Tests.csv")
    df = DataFrame(datatest)

    # Concatenate the data to plot
    datatograph1=hcat(df[:,1],df[:,2],df[:,3],df[:,7]*200)
    datatograph2=hcat(df[:,4],df[:,5],df[:,6],df[:,8])
    #First graph showing the calculation time required to construct the LBLT factorization with the different methods according to the dimension of the problem considered
    plot(df[:,9],datatograph1, label = ["lbl bkaufmann" "lbl bparlett" "lbl rook" "bunchkaufman*200)" ] , ylabel="Temps(s)",xlabel="Dimension", legend=:topleft) #title="Temps de calcul selon la dimension \n (Factorisation)"
    #Second graph showing the calculation time required to solve Ax=b using the previously constructed LBLT factorization according to the dimension of the problem considered
    plot(df[:,9],datatograph2, label = ["lbl bkaufmann solver" "lbl bparlett solver" "lbl rook solver" "Bunchkaufman solver"] , ylabel="Temps(s)",xlabel="Dimension", legend=:topleft) #, title="Temps de calcul selon la dimension \n (Résolution)"
