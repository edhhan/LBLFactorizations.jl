include("pivot_strategies/bkaufmann.jl")
include("pivot_strategies/bparlett.jl")
include("pivot_strategies/rook.jl")


using LinearAlgebra
using Test

#Test for bkaufmann, 4 cases 
@testset begin
    #First case: maximum magnitude of any subdiagonal entry in column 1 = 0, should return pivot_size=1, pivot=(1,1)

    A=[10 0  0  0;
    0  3  0  0;
    0  4  17 0;
    0  1  21 298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot_size-1==0

    #Second case: |a_{11}| >= α * ω_1, should return pivot_size=1, pivot=(1,1) ω_1 is maximum magnitude of any subdiagonal entry in column 1.

    A=[10    0  0  0  ;
       15.6  15 0  0  ;
       1     4  17 0  ;
       0     1  21 298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot_size-1==0

    #Second case: |a_{11}| >= α * ω_1, should return pivot_size=1, pivot=(1,1) ω_1 is maximum magnitude of any subdiagonal entry in column 1.

    A=[10    0  0  0  ;
       15.6  15 0  0  ;
       1     4  17 0  ;
       0     1  21 298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot_size-1==0

    #Third case: |a_{11}|ω_r ≥ α * ω_1^2 , should return pivot_size=1, pivot=(1,1) ω_1 is maximum magnitude of any subdiagonal entry in column 1.

    A=[10    0  0  0  ;
       15.6  15 0  0  ;
       1     4  17 0  ;
       0     1  21 298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot_size-1==0

end
