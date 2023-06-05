restart
load "ZZdFactorizations.m2" 

--shift function
shift = method()
shift(ZZdFactorization) := X -> (
    ZZdfactorization{-(dd^X)_2, -(dd^X)_1}
    )

--direct sum function
directSumMF = method()
directSumMF(ZZdFactorization, ZZdFactorization) := (X, Y) -> (
    --check they are factorizations of same polynomial
    x:= dd^X;
    y:= dd^Y;
    ZZdfactorization{x_1 ++ y_1, x_2 ++ y_2}
    )

--isMatrixFactorization
--returns true if matrices in factorization multiply to f*id in all (2) cyclic permutations
isMatrixFactorization = method()
isMatrixFactorization(ZZdFactorization) := X -> (
    P:= (dd^X)_1*(dd^X)_2; --product of matrices
    S:= (dd^(shift X))_1*(dd^(shift X))_2; --product of shifted matrices
    P == P_(0,0)*id_(source (dd^X)_1) and S == S_(0,0)*id_(source(dd^(shift X))_1)
    )

--EXAMPLES
Q = QQ[x,y]
A = map(Q^2, Q^2, matrix{{x,y}, {-y,x}})
B = map(Q^2, Q^2, matrix{{x,-y}, {y, x}})
X = ZZdfactorization{A,B}

--shift demo
{dd^X, dd^(shift X)}


--boring factorization
i = map(Q^1, Q^1, id_(Q^1))
f = map(Q^1, Q^1, matrix{{x^2 + y^2}})
Y = ZZdfactorization{i, f}

Y' = ZZdfactorization{ map(Q^2, Q^2, matrix{{x, 3*x}, {-y, x*y}}), map(Q^2, Q^2, matrix{{x*y, -3*x}, {y, x}})}

directSumMF(Y, Y')
isMatrixFactorization oo

isMatrixFactorization X
