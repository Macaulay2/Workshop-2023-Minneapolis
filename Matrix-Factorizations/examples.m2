--restart
--needs "ZZdFactorizations.m2"
needs "shift.m2"
--needs "tensorMF.m2"

-------------------------------------
--EXAMPLES
Q = QQ[x,y]
Y = ZZdfactorization{ matrix{{x}}, matrix{{y}}}
dd^oo

A = map(Q^2, Q^2, matrix{{x,y}, {-y,x}})
B = map(Q^2, Q^2, matrix{{x,-y}, {y, x}})
X = ZZdfactorization{A,B}

--shift demo
{dd^X, dd^(shift X)}


--boring factorization
Q = QQ[x,y]
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

dual X
X
D = ZZdfactorization{-dual((dd^X)_2), dual((dd^X)_1)}
dual X

dd^oo

options dual
dualMF = ZZdfactorization{-dual d_1, dual d_2}
maps = {dd^(dualMF)_1, dd^(dualMF)_2}
maps_0*maps_1

homMF = tensorMF(dualMF X, Y)

--mapping cone
f = map(X, X, {
	map(X_0, X_0, m),
	map(X_1, X_1, m),
	map(X_2, X_2, m)}) 

C = cone f
isMatrixFactorization C
dd^C


S = QQ[x,y]
I = ideal(x,y)
x+y
isSubset(ideal(x+y), I)

if not isSubset(ideal(f), I)

instance(x+y, RingElement)
