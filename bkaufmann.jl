using LinearAlgebra
include("max_subdiagonal.jl")

function bkaufmann(A::Hermitian{T}) where T

    α=(1+sqrt(17))/8

    n=size(A,1)

    #μ_1=max |a_ij|

    ω_1=0
    ω_r=0

    pivot=tuple()
    pivot_size=0

    #Looks for the largest absolute value in the subdiagonal element of column 1
    ω_1, r = max_subdiagonal(A, 1)

        a_11=abs(A[1,1])
        
        if (ω_1 == 0)
            pivot = (0,0)
            pivot_size = 0

        elseif a_11 >= α*ω_1
  
            pivot_size=1
            #The pivot is a_11 is first element of the tuple pivot
            pivot=(1,0)
            
        elseif (a_11*(ω_r = max_subdiagonal(A, r)[1])  >= α*ω_1^2) 

            pivot_size=1
            #The pivot is a_11 is first element of the tuple pivot
            pivot=(1,0)

        elseif (a_rr = abs(A[r,r])) >= α*ω_r

            pivot_size=1
            #The pivot is a_rr is first element of the tuple pivot
            pivot=(r,0)

        else
            
            pivot_size=2
            pivot=(1,r)

        end

    return pivot, pivot_size
end
