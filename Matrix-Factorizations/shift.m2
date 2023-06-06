restart
load "ZZdFactorizations.m2" 

--shift function
shift = method()
shift(ZZdFactorization) := ZZdFactorization => X -> (
    ZZdfactorization{-(dd^X)_2, -(dd^X)_1}
    )

--direct sum function
directSumMF = method()
directSumMF(ZZdFactorization, ZZdFactorization) := ZZdFactorization => (X, Y) -> (
    --check they are factorizations of same polynomial
    x:= dd^X;
    y:= dd^Y;
    ZZdfactorization{x_1 ++ y_1, x_2 ++ y_2}
    )

--isMatrixFactorization
--returns true if matrices in factorization multiply to f*id in all (2) cyclic permutations
isMatrixFactorization = method()
isMatrixFactorization(ZZdFactorization) := X -> (
    P:= matrix entries ((dd^X)_1*(dd^X)_2); --product of matrices
    S:= matrix entries ((dd^(shift X))_1*(dd^(shift X))_2); --product of shifted matrices
    P == matrix entries (P_(0,0)*id_(source (dd^X)_1)) and S == matrix entries (S_(0,0)*id_(source(dd^(shift X))_1))
    )

--function for generating matrix factorization from a
--module M over hypersurface ring
tailMF = method()
tailMF(Module) := ZZdFactorization => M -> (
    F := res(M, LengthLimit=>dim ring M + 1);
    C := chainComplex sub(F.dd_(dim ring M + 1), ambient ring M);
    h := (nullhomotopy extend(C, C, ((ring M).relations)_(0,0)*id_(C_0)))_0;
    ZZdfactorization{sub(F.dd_(dim ring M + 1), ambient ring M), h}
    )

--cone
--in: map of MFs
--out: new MF
cone(ZZdFactorizationMap) := ZZdFactorization => f -> (
    s := source f;
    t := target f;
    m1 := matrix{{(dd^s)_1, map(s_0, t_0, 0)},
    {f_1, -(dd^t)_2}};
    m2 := matrix{{(dd^s)_2, map(s_1, t_1, 0)}, 
    {f_0, -(dd^t)_1}};
    ZZdfactorization{m1, m2}
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

--tailMF demo
R = QQ[x,y]/ideal(x^3 + y^3)
M = comodule ideal vars R
X = tailMF M
S = ring X
m = (x+y)*id_(S^2)

--mapping cone
f = map(X, X, {
	map(X_0, X_0, m),
	map(X_1, X_1, m),
	map(X_2, X_2, m)}) 

C = cone f
isMatrixFactorization C
