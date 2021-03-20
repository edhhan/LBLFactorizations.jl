using LinearAlgebra
include("max_subdiagonal.jl")

"""
Rook pivoting strategy
"""
function rook(A::Hermitian{T}) where T

         
    α = (1+sqrt(17))/8
    n = size(A,1)

    #μ_1=max |a_ij|

    ω_i = 0
    ω_r = 0

    pivot = tuple(0,0)
    pivot_size = 0

    #Looks for the largest absolute value in the subdiagonal element of column 1
    ω_i, r = max_subdiagonal(A, 1)

        a_11 = abs(A[1,1])
        
        if (ω_i == 0)
            pivot = (0,0)
            pivot_size = 0

        elseif a_11 >= α*ω_i
  
            pivot_size = 1
            #The pivot is a_11 is first element of the tuple pivot
            pivot = (1,0)
            
        else
            i = 1
            while(pivot_size==0)

               ω_r, r = max_subdiagonal(A, i)
                #Voir ici ce qu'on fait si on arrive à la fin sans avoir trouver de pivot, est-ce qu'on prend le dernier en mémoire ?
               if r == 0
                    pivot = (i,0)
                    pivot_size = 1

               elseif abs(A[r,r]) >= α*ω_r 

                pivot = (r,0)
                pivot_size = 1

               elseif ω_i == ω_r

                pivot= (i,r)
                pivot_size = 2

               else

                i=r
                ω_i = ω_r

               end

            end

        end

    return pivot, pivot_size

end