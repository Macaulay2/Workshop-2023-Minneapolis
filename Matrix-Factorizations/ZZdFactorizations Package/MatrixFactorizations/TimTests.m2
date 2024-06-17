restart
needs "ZZdFactorizations.m2"
needs "MF_functions.m2"


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

S = QQ[x]
A = matrix{{1,1}};
B = matrix{{0},{0}};
X = ZZdfactorization{A,B}
isdFactorization(X)

///

TEST ///

S = QQ[x,y]
X = ZZdfactorization{x,x}
Y = ZZdfactorization{y,y}
Z = X++Y
assert((isdFactorization(Z))_0 == false)

///

TEST ///

S = QQ[x,y,z]
X = ZZdfactorization{S^3,S^5,S^2,S^1}
isdFactorization X
///