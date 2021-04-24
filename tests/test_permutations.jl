using LinearAlgebra
using Test

n=5
A=rand(Float64, n,n).*100
AL = Hermitian(A,:L)
AU = Hermitian(A,:U)
display(AL)
display(AU)
AL.data[3,3]=AL.data[4,4]
display(AL)
display(AU)
print(AU.uplo)
#=
#s=premier membre de pivot,p=deuxieme membre de pivot
if (AU.uplo==:U)

    temp=A.data[s,s+1:p-1]
    A.data[s,s+1:p-1]=A.data[s+1:p-1,p]'
    A.data[s+1:p-1,p]=temp'

    A.data[s,p]=A.data[s,p]'

    temp=A.data[s,p+1:end]
    A.data[s,p+1:end]=A.data[p,p+1:end]
    A.data[p,p+1:end]=temp

    temp=A.data[1:s-1,s]
    A.data[1:s-1,s]=A.data[1:s-1,p]
    A.data[1:s-1,p]=temp   


    temp=A.data[s,s]
    A.data[s,s]=A.data[p,p]
    A.data[p,p]=temp   
end

if (A.uplo==:L)

#Même logique qu'en haut mais inversé
end
=#
#=
# Permutation inplace on lines 
temp = A_prime[p[1], :]
A_prime[p[1], :] = A_prime[p[2], :]
A_prime[p[2], :] = temp
                   
# Permutation inplace on columns
temp = A_prime[:, p[1]]
A_prime[:, p[1]] = A_prime[:,p[2]]
A_prime[:, p[2]] = temp
=#