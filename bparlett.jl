using LinearAlgebra
include("max_subdiagonal.jl")

function bparlett(A::Hermitian{T}) where T

    α=(1+sqrt(17))/8

    n=size(A,1)

    #μ_1=max |a_ij|

    μ_0=0
    μ_1=0

    μ_1_r=0

    μ_0_p=0
    μ_0_q=0

    pivot=tuple()
    pivot_size=0

    for j in 1:n
        for i in j:n

            v=abs(A[i,j])

            #Max dans la diagonale
            if (i==j)
                if(v>μ_1)
                    μ_1=A[i,j]
                    μ_1_r=i
                end
            end
            #Max de tous les éléments
            if(v>μ_0) 
                μ_0=A[i,j]
                μ_0_p=i
                μ_0_q=j
            end
        end
        
        if (μ_1 >= α*μ_0)
            pivot_size=1
            pivot=(μ_1_r,0)
        else
            pivot_size=2
            pivot=(μ_0_p,μ_0_q)            
        end



    end

    return pivot, pivot_size
end