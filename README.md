# julia-factorization-bunch-kaufman

# LBL: Bloc Diagonal Factorization for Hermitian Matrices

IMPORTANT: This package isn't efficient compared to bunchkaufman() and was done for academic purpose only. Use bunchkaufman() from LinearAlgebra for optimal results.

An implementation of Bloc diagonal factorization

Please cite this repository if you use lbl.jl in your work

This package is appropriate for indefinite hermitian matrix A where we want a factorization of the form PAP*=LBL*.

lbl.jl should not be expected to be as fast as bunchkaufman.jl which uses LAPACK.

The main advantages of lbl.jl are that it is implemented in Julia so it does not require external compiled dependencies and it can work with all the julia data types.

# Installing

```julia
julia> ]
pkg> add lbl
```

# Usage

The only exported functions are `lbl`, `lbl_solve`.


`lbl` returns a factorization in the form of a LBL object.
The `lbl_solve` method is implemented for objects of type `LBL` so that
solving a linear system is as easy as
```julia
LBLT = lbl(A,strategy)  # LBLᵀ factorization of A with pivoting strategy strategy: "rook", "bkaufmann" or "bparlett"

x = lbl_solve(LBLT, b) # solves Ax = b
```

Factors can be accessed as `LBLT.L`, `LBLT.B`, and the permutation vector as `LBLT.permutation`.
L factor is unit lower triangular and the block diagonal matrix B is tridiagonal and contains independant submatrix 2x2 and 1x1.

# References

The strategies and the factorization were based on:

N. J. Higham,Accuracy and stability of numerical algorithms.  SIAM, 2002
