needs "ZZdFactorizations.m2"
needsPackage "TensorComplexes" --Koszul factorization method needs this 

--shift function
--Only works in d = 2 case
shift = method()
shift(ZZdFactorization) := ZZdFactorization => X -> (
    ZZdfactorization{-(dd^X)_2, -(dd^X)_1}
    )

--direct sum function
--Only works in d = 2 case
directSumMF = method()
directSumMF(ZZdFactorization, ZZdFactorization) := ZZdFactorization => (X, Y) -> (
    --check they are factorizations of same polynomial
    x:= dd^X;
    y:= dd^Y;
    ZZdfactorization{x_1 ++ y_1, x_2 ++ y_2}
    )


--In: map of MFs
--Out: new MF
--Only works in d=2 case
cone(ZZdFactorizationMap) := ZZdFactorization => f -> (
    s := source f;
    t := target f;
    m1 := matrix{{(dd^s)_1, map(s_0, t_0, 0)},
    {f_1, -(dd^t)_2}};
    m2 := matrix{{(dd^s)_2, map(s_1, t_1, 0)}, 
    {f_0, -(dd^t)_1}};
    ZZdfactorization{m1, m2}
    )


--Checks if differential matrices (and all cyclic shifts of them) multiply to same polynomial
--If so, it will output that polynomial
isdFactorization = method()
isdFactorization(ZZdFactorization) := X -> (
    P:= matrix entries product(apply(X.period, i->(dd^X)_i));
    f:= P_(0,0);
    S:= for i to X.period-1 list(
	product(apply(X.period, j->(dd^X)_(i+j)))
	); --S is a list of product of all cyclic shifts of matrix multiplication
    --S)
    if same S then (same S, f) else (same S, "no potential")
    --error "Your ZZdFactorization is not an actual matrix factorization."
    )

-- function for generating matrix factorization from a
-- module M over hypersurface ring
-- output is a d=2 factorization
tailMF = method()
tailMF(Module) := ZZdFactorization => M -> (
    F := res(M, LengthLimit=>dim ring M + 1);
    C := chainComplex sub(F.dd_(dim ring M + 1), ambient ring M);
    h := (nullhomotopy extend(C, C, ((ring M).relations)_(0,0)*id_(C_0)))_0;
    ZZdfactorization{sub(F.dd_(dim ring M + 1), ambient ring M), h}
    )

--note: dual is a method with options
--d = 2 version
dual(ZZdFactorization) := ZZdFactorization => {} >> o -> X -> (
    if not(X.period == 2) then error "since period of X is not 2, please input a root of unity";
    ZZdfactorization{-dual((dd^X)_2), dual((dd^X)_1)}
    )
-- d>2 requires one to input a dth root of unity
dual(ZZdFactorization, RingElement) := ZZdFactorization => {} >> o -> (X,omega) -> (
    ZZdfactorization for i to X.period-1 list (-omega^i)*dual((dd^X)_(X.period - i))
    )

--d=2 version
Hom(ZZdFactorization, ZZdFactorization) := ZZdFactorization => (X,Y) -> (
    tensorMF(dual X, Y)
    )

--d>2 requires one to input a dth root of unity
Hom(ZZdFactorization, ZZdFactorization, RingElement) := ZZdFactorization => (X,Y,omega) -> (
    dTensor(dual(X,omega),Y,omega)
    )

--d=2 case
tensorMF = method()
tensorMF(ZZdFactorization,ZZdFactorization) := (X,Y) -> (
    if not(X.period==2 and Y.period==2) then error "Both inputs should have period 2. Use dTensor method instead.";
    A1 := id_(X_0) ** (dd^Y)_1;
    B1 := (dd^X)_1 ** id_(Y_0);	   
    C1 := (dd^X)_2 ** id_(Y_1);
    D1 := -id_(X_1)** (dd^Y)_2;
    diff1 := matrix{{A1, B1},{C1,D1}};
    A2  := id_(X_0) ** (dd^Y)_2;
    B2 := (dd^X)_1 ** id_(Y_1);
    C2 := (dd^X)_2 ** id_(Y_0);
    D2 := -id_(X_1)** (dd^Y)_1;
    diff2 := matrix{{A2, B2},{C2,D2}};
    ZZdfactorization{diff1, diff2}
)

--d>2 requires one to input a root of unity
dTensor = method()
dTensor(ZZdFactorization, ZZdFactorization,RingElement) := (F,G,t) -> (
    --put in check for F.period = G.period
    N := F.period;
    dF := dd^F;
    dG := dd^G;
    M := for k to N-1 list(
	for i to N-1 list(
	    for j to N-1 list(
		if i == j then (t^i*id_(F_i))**(dG_(k-i+1))		
		else if (j == (i+1)%N) then (dF_(i+1))**(id_(G_((k-i)%N)))
		else 0
	    )));
    ZZdfactorization(for i to #M-1 list matrix M_i)
)

-----KoszulMF
--koszul matrix factorization
--in: list of even number of elements in ring
--out: matrix factorization of sum products of pairs (L#i)*(L#(i+1))

koszulMF = method()
koszulMF(List) := L -> (
    X := ZZdfactorization{L#0, L#1};
    for i from 1 to (#L)//2-1 do X = tensorMF(X, ZZdfactorization{ L#(2*i), L#(2*i+1)});
    X
)
-- 2 case is not consistent breaks a list into 2, rather than inputs 2 lists


---inputs a list of d lists each of length n outputs the tensor product of n ZZdfactorizations
koszulMF(List, RingElement) := (L, omega) -> (
    n := #(L#0);
    d := #L;
    Z := for i to n-1 list ZZdfactorization for j to d-1 list (L#j)#i;
    T := Z#0;
    for i from 1 to n-1 do T = dTensor(T,Z#i, omega);
    T
)


------KoszulMF from ideal and function
koszulMFf = method()
koszulMFf(List, RingElement) := (L,f) -> (
    I := ideal L;
    M := f //gens(I);
    N := transpose (gens I) | M;     
    E  := entries N;
    F := select(E, i -> not(i#1 == 0)); 
    A := flatten F;
    koszulMF(A)
)

koszulMFf(Ideal, RingElement) := (I,f) -> (
    koszulMFf(flatten entries gens I, f)
)

eulerMF = f -> (koszulMFf(ideal jacobian f , f))

--In: List L generating an ideal I such that polynomial f lies in I^{d-1}, where d is the desired period of factorization
--In: root of unity omega and integer e=d-1 such that omega^(e+1)=1
--Out: factorization of f with period d
--omega should be an (e+1) root of unity!
koszulMFf(List, RingElement, ZZ, RingElement) := (L,f, e, omega) -> (
    if not(omega^(e+1) == sub(1,ring omega)) then error "ring element must be a d^th root of unity where d is the third input";
    G := matrix {L};
    Gd := multiSubsets(flatten entries G, e);
    prods := for i to #Gd-1 list product(Gd_i);
    M := f //matrix{prods};
    N := for i to #Gd-1 list flatten {Gd#i, (entries M)#i};
    F := select(N, i -> not(last i == 0));
    koszulMF(entries transpose matrix F, omega)
	)

--i = homological degree
--j = F_j degree
freeMods = method()
freeMods(ZZdFactorization, ZZdFactorization) := (F,G) -> (
    N:= F.period;
for i to N-1 list(
for j to N-1 list(
    {j, (N + i - j)%N} => (F_j)**(G_((N+i-j)%N))
    )))

adjoinRoot = (d,Q) -> (S1 := Q[t];
    factored := factor(t^(d)-1);
    cyclo := (factored#(#factored-1))#0;
    S := S1/(cyclo)
    )

--this gives a trivial d-fold factorization of a monomial
trivialFactorization = method()
trivialFactorization(RingElement) := f -> (if not(#(terms f)==1) then error "Expected ring element to be monomial";
    theDiffs = (flatten((toList factor(f))/(i->toList(i#1:i#0))))/(j->matrix{{j}});
    ZZdfactorization theDiffs
    )
