
TEST ///

assert(1==1)

///


TEST ///

--zero factorization on free modules of rank 1
S = ZZ/101[x];
assert(ZZdfactorization{0_S, 0_S} == ZZdfactorization{matrix{{0_S}}, matrix{{0_S}}})
assert(ZZdfactorization{0_S, 0_S} == ZZdfactorization{S^1, S^1})

--direct sum works as expected
I = ZZdfactorization{1_S, 1_S};
J = ZZdfactorization{x, 2_S};
assert(I++I == ZZdfactorization{id_(S^2), id_(S^2)})
X = I++J;

assert(prune image X_[0] == I)
assert(prune image X_[1] == J)
assert(prune image X^[0] == I)
assert(prune image X^[1] == J)
assert(prune coker X_[0] == J)
assert(prune coker X_[1] == I)

///


TEST ///

S = ZZ/101[x,y,z];
X = ZZdfactorization{x, y, z};
--identity and 0 map are well defined, commutative
assert isWellDefined id_X
assert isCommutative id_X
assert isWellDefined(0*id_X) --error: key not found in hash table
assert isCommutative(0*id_X)
--isdFactorization works
assert((isdFactorization X)_0)
--diff of End composes to 0
assert( ((End adjoinRoot(X,t)).dd)^3 == 0)
--isdFactorization returns true and factorization of 0
assert (isdFactorization(End adjoinRoot(X, t)))_0
assert((isdFactorization(End adjoinRoot(X, t)))_1 == 0)
--isZZdComplex
assert isZZdComplex End adjoinRoot(X, t)

    
Y = ZZdfactorization{x,y};
--End of length 2 factorization is a complex
assert(((End(Y)).dd)^2 == 0)
assert (isdFactorization End Y)_0
assert ((isdFactorization End Y)_1 == 0)
assert isZZdComplex End Y

///

TEST ///

--shift composes well
S = ZZ/101[x,y,z];
X = ZZdfactorization{x,y,z};
X' = adjoinRoot(X, t);
assert((X'[1])[1] == X'[2])
assert( X'[3] == X')

--periodicity of modules and differential
assert(X_0 == X_18)
assert(X.dd_1 == X.dd_31)

--checking that id and 0 work well
assert( id_X == 1)
assert(0*id_X == 0)
assert(id_X - id_X == 0)
assert(isCommutative(id_X - id_X))

///


TEST ///

--folding complex gives factorization of zero
S = ZZ/101[x,y,z];
K = koszulComplex vars S;
F = Fold(K, 2);
assert( (F.dd)^2 == 0)

///

TEST ///

--tailMF gives a d factorization
R = ZZ/101[x,y,z]/ideal(x^2 + y^2 + z^2);
M = coker vars R;
assert( (isdFactorization(tailMF M))_0)
assert( (isdFactorization(tailMF M))_1 == (ideal R)_0)

///

TEST ///

S = QQ[x]
A = matrix{{1,1}};
B = matrix{{0},{0}};
X = ZZdfactorization{A,B}
isdFactorization(X)

///

TEST ///

S=QQ[x]
X = ZZdfactorization{x,x,x}
a = map(X,X,id_X)
assert(a_0 == id_(S^1))
assert(a_1 == id_(S^1))
assert(source a_0 == source a_2)

///

TEST ///

S = ZZ/13[x,y]
X = ZZdfactorization{x,y}
assert((isdFactorization(X))_0)
assert(X.period%2 == 0)
assert(X.dd_0*X.dd_1==X.dd_1*X.dd_2)
assert(X.dd_0 == X.dd_2)
assert(X.dd_1 == X.dd_(-1))
assert(X_0 == X_2)
assert(X_1 == X_3)

///

TEST ///

S = ZZ/7[x,y,z]
X = ZZdfactorization{x,y,z}
assert((isdFactorization(X))_0)
assert(X.period%3 == 0)
assert(X.dd_0*X.dd_1*X.dd_2==X.dd_1*X.dd_2*X.dd_3)
assert(X.dd_0 == X.dd_3)
assert(X.dd_1 == X.dd_(-2))
assert(X_0 == X_3)
assert(X_1 == X_4)

///

TEST ///

S = ZZ/19[x,y,t]/(t^2+t+1)
X = ZZdfactorization{x,x,x}
Y = ZZdfactorization{y,y,y}
Z = dTensor(X,Y,t)
assert(isdFactorization(dTensor(X,X,t))_1==2*x^3)
assert((isdFactorization(Z))_1 == x^3+y^3)
assert(t^3==1)

///

TEST ///

S = QQ[x]
A = matrix{{1,1}};
B = matrix{{0},{0}};
X = ZZdfactorization{A,B}
isdFactorization(X)

///

TEST ///

--some basic tests on maps
S = ZZ/101[x,y,z];
X = ZZdfactorization{x,y,z};

assert((id_X)_1 == (id_X)_31)
assert( (id_X)^2 == (id_X))

///

