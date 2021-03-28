include("../pivot_strategies/bparlett.jl")
include("../pivot_strategies/bkaufmann.jl")
include("../pivot_strategies/rook.jl")


using LinearAlgebra
using Test

#Test for bkaufmann, 5 cases 
@testset begin
    #First case: maximum magnitude of any subdiagonal entry in column 1 = 0, should return pivot_size=0, pivot=-1

    A=[-10 0  0    0  ;
        0  15  0   0  ;
        0  -4 -17  0  ;
        0   1  33 -298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1]==-1 && pivot_size==0

    #Second case: |a_{11}| >= α * ω_1, should return pivot_size=1, pivot=(1,1) ω_1 is maximum magnitude of any subdiagonal entry in column 1.

    A=[-10    0  0   0  ;
        15.6  15 0   0  ;
        1    -4 -17  0  ;
        0     1  33 -298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot_size==1

    #Third case: |a_{11}| * ω_r >= α * ω_1^2, should return pivot_size=1, pivot=(1,1) ω_r is maximum magnitude of any OFF diagonal entry in column r.
    # r is row index of first (SUBdiagonal) entry of maximum magnitude in column 1. In this case: |a_{11}|=-10, ω_r=|-33|, ω_1=20. Should return pivot=(1,1)

    A=[-10   0   0    0  ;
       15.6  15  0    0  ;
       20   -33  -17   0  ;
       2     1   32  -298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot_size==1

    #Fourth case: |a_{rr}| ≥ α * ω_r , should return pivot_size=1, pivot=(1,r)=(1,3) ω_1 is maximum magnitude of any subdiagonal entry in column 1.

    A=[-10   0   0    0  ;
       15.6  15  0    0  ;
       20   -2  -17   0  ;
       2     1   2  -298]

    A=Hermitian(A,:L)

    pivot, pivot_size = bkaufmann(A)

    @test pivot[1][1]==1 && pivot[1][2]==3 && pivot_size==1

    #Fifth case: if all of the above didn't work , should return pivot_size=2, pivot=(2,r) ω_1 is maximum magnitude of any subdiagonal entry in column 1.

    A=[-0.01   0   0    0  ;
       15.6  15  0    0  ;
       20   -15  -17   0  ;
       2     1   30  -298]

    A=Hermitian(A,:L)
    
    pivot, pivot_size = bkaufmann(A)
 
    @test pivot[1][1]==2 && pivot[1][2]==3  && pivot_size==2   


#Test for bparlett, 2 cases


    #First case: μ_1 >= α * μ_0 where μ_0 is maximum value of matrix, μ_1 is maximum value of diagonal. Should return pivot_size=1, pivot=(1,r). r is index of diagonal value. In this case, r=3

    A=[ 1   0   0     0  ;
        0   0   0     0  ;
        465  -4 -298  0  ;
        298   1  33   297]

    A=Hermitian(A,:L)

    pivot, pivot_size = bparlett(A)

    @test pivot[1][1]==1 && pivot[1][2]==3 && pivot_size==1

    #Second case: if above is not true: μ_1 < α * μ_0 where μ_0 is maximum value of matrix, μ_1 is maximum value of diagonal. Should return pivot_size=2, pivot=(1,p) and (2,q). p is row index of maximum value in upper triangular. q is column index in upper triangular. In this case, p=3, q=4

    A=[ 1   0   0      0  ;
        0   0   0      0  ;
        465  -4 -298   0  ;
        298   1  466   297]

    A=Hermitian(A,:L)

    pivot, pivot_size = bparlett(A)
    println(pivot, pivot_size)
    @test pivot[1][1]==1 && pivot[1][2]==3 && pivot[2][1]==2 && pivot[2][2]==4 && pivot_size==2



#Test for rook, for 4x4 matrices, 3 cases

    #First case: ω_1=0 , should return pivot_size=0, pivot=-1

    A=[ -10  0   0     0  ;
         0   0   0     0  ;
         0  -4  -298   0  ;
         0   1   33    297]

    A=Hermitian(A,:L)

    pivot, pivot_size = rook(A)

    @test pivot[1]==-1 && pivot_size==0

    #Second case: |a_{11}| >= α * ω_1, should return pivot_size=1, pivot=(1,1)

    A=[ -10     0   0     0  ;
         1      0   0     0  ;
        -15.6  -4  -298   0  ;
         15.4   1   33    297]

    A=Hermitian(A,:L)

    pivot, pivot_size = rook(A)

    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot_size==1
    
    #Third case: First iteration first case: |a_{rr}| >= α * ω_r, i=1, r=3, should return pivot_size=1, pivot=(1,r)=(1,3)

    A=[ 1   0   0      0  ;
        1   0   0      0  ;
        -20  -4 -298   0  ;
        15.4   1  40   297]

    A=Hermitian(A,:L)

    pivot, pivot_size = rook(A)
    println(pivot, pivot_size)
    @test pivot[1][1]==1 && pivot[1][2]==3 && pivot_size==1

    #Third case: First iteration second case: ω_i = ω_r, i=1, r=3, should return pivot_size=2, pivot=(1,i)=(1,1) and (2,r)=(2,3)

    A=[  1      0   0    0  ;
         1      0   0    0  ;
        -20   -4   -10   0  ;
         15.4   1   20   297]

    A=Hermitian(A,:L)

    pivot, pivot_size = rook(A)
    println(pivot, pivot_size)
    @test pivot[1][1]==1 && pivot[1][2]==1 && pivot[2][1]==2 && pivot[2][2]==3 && pivot_size==2

    #Third case: Second iteration first case: |a_{rr}| >= α * ω_r, i=3, r=4, ω_r= 0, |a_{rr}|=297 should return pivot_size=1, pivot=(1,r)=(1,4)

    A=[  1      0   0    0  ;
         1      0   0    0  ;
        -20   -4   -10   0  ;
         15.4   1   30   297]

    A=Hermitian(A,:L)

    pivot, pivot_size = rook(A)
    println(pivot, pivot_size)
    @test pivot[1][1]==1 && pivot[1][2]==4  && pivot_size==1
    
 

end