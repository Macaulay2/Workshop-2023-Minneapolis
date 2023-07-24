restart
needs "ZZdFactorizations.m2"
needs "functionsMF_new.m2"
Q = QQ[x,y,z]
R = Q/(x^3 + y^3+z^3)
M = R^1/ideal(x^2,y,z)

tailMF = M -> (
    R:=ring M;
    F := res(M,LengthLimit=>dim R +1);
    C := chainComplex(sub(F.dd_(dim R +1), ambient R));
    h := (nullhomotopy extend(C,C,(R.relations)_(0,0)*id_(C_0)))_0;
    X := ZZdfactorization{sub(F.dd_(dim R+1),ambient R),h}	
    )

X = tailMF M

mcmRank = X -> (
    R := (ring X)/ideal(equationMF X);
    use R;
    (rank (R**coker (dd^X)_1),rank (R**coker (dd^X)_2))
    )

mcmRank X

isUlrich = X -> (
   (rank X_0 == ((rankMF X)_0)*(first degree equationMF X),rank X_1 == ((rankMF X)_1)*(first degree equationMF X))
    )	 

isUlrich X

-- input a dMF X and an integer k
-- output the kth shift of X
-- can we get it to automatically shift by 1 if you don't input an integer? 
shift = method()
shift(ZZdFactorization,ZZ) := (X,k) -> (
    Y := for i to X.period-1 list(
	(dd^X)_(i+k+1)
	);
    ZZdfactorization(Y)
    )

S = QQ[x,y,z,w]
Z = ZZdfactorization{x,y,z,w}
s = shift(Z,2)
dd^s

-- direct sum of two d-fold MFs of the same polynomial
directSumMF = method()
directSumMF(ZZdFactorization,ZZdFactorization) := (X,Y) -> (
    -- add checks (X,Y should be MFs of same period and of same polynomial)
    x := dd^X;
    y := dd^Y;
    z := for i from 1 to X.period list(
	x_i ++ y_i
	);
    ZZdfactorization(z)
    )

m = directSumMF(Z,shift(Z,1))
isdFactorization(m)
dd^m


-- input: dMF and an integer k
-- output: the (d-1)MF which has d_k*d_{k+1} now in the kth spot 
collapseMF = method()
collapseMF(ZZdFactorization , ZZ) := (X,k) -> (
    if X.period == 2 then error;
    s := shift(X,k-1);
    L := for i from 3 to X.period list(
	(dd^s)_i
	);
    x := ZZdfactorization({(dd^s)_1*(dd^s)_2}|L);
    shift(x, X.period - k)
    )

S = QQ[x_1..x_6]
Z = ZZdfactorization{x_1,x_2,x_3,x_4,x_5,x_6}
dd^Z
C = collapseMF(Z,4)
dd^C

-- input: dMF and integers n and r
-- output: 2MF obtained by composing r maps starting at d_n
-- (second matrix is the composition of the rest of the maps)
fullCollapse = method()
fullCollapse(ZZdFactorization,ZZ,ZZ) := (X,n,r) ->(
    	if X.period == 2 then error;
	x := product(apply(r, i->(dd^X)_(n+i)));
	y := product(apply(X.period-r, i->(dd^X)_(n+i+r)));
	ZZdfactorization{x,y}
    )

Y = fullCollapse(Z,3,3)
dd^Y



-- input: dMF X (in theory a dMF of 0) and integers n and r
-- output: the complex obtained by composing r maps of X starting at degree n
toComplex = method()
toComplex(ZZdFactorization,ZZ,ZZ) := (X,n,r) -> (
    -- add same checks that David did
    c := fullCollapse(X,n,r);
    L := for i to 2 list((dd^c)_(i+1));
    chainComplex(L)
    )

c = toComplex(Z,1,3)
c.dd

-- compute rth layer of homology at homological degree n (or n \pm 1?)
-- these homologies will always be zero if input a dMF which has been tensored to the hypersurface (see x^5 example below)
-- they will be possibly non-zero if you input a random dMF of 0 (for example tensor a factorization of f with one of -f)
homologyMF = method()
homologyMF(ZZdFactorization,ZZ,ZZ) := (X,n,r) -> (
    C := toComplex(X,n,r);
    h := HH C;
    {h_1,h_2}
    )

T = QQ[x]/(x^5)
X = ZZdfactorization{x,x,x,x,x}
h = homologyMF(X,1,4)
h_1 == 0


-- input: dMF X and integers n and r
-- output: list of lengths of homologies of "depth" r starting at homological position n.
-- The n may not really be necessary here.
-- What I want to do here, is add up this list as a weighted sum a + wb + w^2c + ... where w is a 
-- primitive X.period root of unit.
ecMF = method()
ecMF(ZZdFactorization,ZZ,ZZ) := (X,n,r) -> (
    for i to X.period-1 list(
	degree (homologyMF(X,n+i,r))_0
	)
    )

S = QQ[x,y,t]/(t^2+t+1)
X = ZZdfactorization{x+y,x+t*y,x+t^2*y}
Y = ZZdfactorization{-(x+y),-(x+t*y),-(x+t^2*y)}
isdFactorization(X)
Z = dTensor(X,Y,t)
dd^Z
ecMF(Z,1,1)
ecMF(Z,1,2)

-- Example: f=xyz
-- Output list of ecMF doesn't always have to be equal! (guess: for isolated singularites,
-- they are usually equal, i.e, the weigted sum is usually zero)
S = QQ[x,y,z,t]/(t^2+t+1)
X = ZZdfactorization{x,y,z}
Y = ZZdfactorization{-y,-x,-z}
Z = dTensor(X,Y,t)
dd^Z
ecMF(Z,1,1)
ecMF(Z,1,2)

S = QQ[x,y,z,u,v,w,t]/(t^2+t+1)
X1 = ZZdfactorization{x,y,z}
X2 = ZZdfactorization{u,v,w}
Z1 = dTensor(X1,X2,t)
Y1 = ZZdfactorization{-x,-z,-y}
Y2 = ZZdfactorization{-w,-u,-v}
Z2 = dTensor(Y1,Y2,t)
Z = dTensor(Z1,Z2,t)
isdFactorization(Z1)
ecMF(Z,1,1)
ecMF(Z,1,2)


S = QQ[x,y,t]/(t^4+t^3+t^2+t+1)
X1 = ZZdfactorization{x,x,x,x,x}
X2 = ZZdfactorization{y,y,y,y,y}
X = dTensor(X1,X2,t)
Y = ZZdfactorization{-(x+y),-(x+t*y),-(x+t^2*y),-(x+t^3*y),-(x+t^4*y)}
isdFactorization(Y)
Z = dTensor(X,Y,t)
ecMF(Z,1,1)
ecMF(Z,1,2)
--don't try to run ecMF(Z,1,2)...
-- bad example anyway. X is just a direct sum of copies of shifts of Y
