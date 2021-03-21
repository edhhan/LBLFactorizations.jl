using LinearAlgebra

function PermuteMatrix(A,pivot_array)

    n=length(pivot_array)
    n_A=size(A,1)
    s=1

    for i=1:n
        P = Matrix(1.0*I, n_A, n_A)
        k=0
        for p in pivot_array[s]

            if p!=-1
            idx1 = p[1]
            idx2 = p[2]

            # Permutation
            temp = P[:,s+idx1-1]
            P[:, s+idx1-1] = P[:, s+idx2-1]
            P[:, s+idx2-1] = temp
            end
            k+=1
        end
        s+=k
        A=P*A*P'
        
    end

end