include("pivot_strategies/bkaufmann.jl")
include("pivot_strategies/bparlett.jl")
include("pivot_strategies/rook.jl")

using LinearAlgebra


function permutation_matrix(pivot,n)

    P=Matrix(1.0I,n,n)
    P[pivot[1], pivot[1]]=0
    P[pivot[2], pivot[2]]=0

    P[pivot[2], pivot[1]]=1
    P[pivot[1], pivot[2]]=1

    return P
    
end


function LDL_MGD(A)
    n_A=size(A,1)
    L=zeros(n_A,n_A)
    BD=zeros(n_A,n_A)
    Permutations=[]

    s=1
    A_prime=copy(A)
    n=size(A_prime,1)

    while(s<n_A)

        pivot, pivot_size = bkaufmann(A_prime)
        if pivot_size!=0
            for i = 1:pivot_size
                P=permutation_matrix(pivot[i],n)
                A_prime=P*A_prime*P'
                push!(Permutations, (s+pivot[i][1]-1, s+pivot[i][2]-1))
            end
        end
        println(pivot_size)
        #Schur
        if pivot_size!=0
            B=A_prime[(pivot_size+1):end,(pivot_size+1):end]
            C=A_prime[(pivot_size+1):end,1:pivot_size]
            E=A_prime[1:pivot_size,1:pivot_size]

            A_prime=B-C*inv(E)*C'
            n=size(A_prime,1)
            L[s:end,s:s+pivot_size-1]=vcat(Matrix(1.0*I, pivot_size, pivot_size), C*inv(E) )
            BD[s:s+pivot_size-1,s:s+pivot_size-1]=E
            s+=pivot_size


        else 
            pivot_size=1 #Théoriquement, on pourrait le merger en haut en faisant if pivot_size=0, pivot_size=1
            B=A_prime[(pivot_size+1):end,(pivot_size+1):end]
            C=A_prime[(pivot_size+1):end,1:pivot_size]
            E=A_prime[1:pivot_size,1:pivot_size]
            
            A_prime=B-C*inv(E)*C'
            n=size(A_prime,1)
            L[s:end,s]=vcat(Matrix(1.0*I, pivot_size, pivot_size), C*inv(E) )
            BD[s:s+pivot_size-1,s:s+pivot_size-1]=E
            s+=pivot_size       


        end


    end
    if s==n_A
        L[n_A,n_A]=1
        BD[n_A,n_A]=A_prime[1,1]
    end

    return L, BD, Permutations
end



