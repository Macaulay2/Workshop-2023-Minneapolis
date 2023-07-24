needs "ZZdFactorizations.m2"
needsPackage "TensorComplexes" --Koszul factorization method needs this 

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
-*isMatrixFactorization = method()
isMatrixFactorization(ZZdFactorization) := X -> (
    P:= matrix entries ((dd^X)_1*(dd^X)_2); --product of matrices
    S:= matrix entries ((dd^(shift X))_1*(dd^(shift X))_2); --product of shifted matrices
    P == matrix entries (P_(0,0)*id_(source (dd^X)_1)) and S == matrix entries (S_(0,0)*id_(source(dd^(shift X))_1))
    )*-

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

--dual
--note: dual is a method with options
-*dual(ZZdFactorization) := ZZdFactorization => {} >> o -> X -> (
    if not(X.period == 2) then error "since period of X is not 2, please input a root of unity";
    ZZdfactorization{-dual((dd^X)_2), dual((dd^X)_1)}
    )

dual(ZZdFactorization, RingElement) := ZZdFactorization => {} >> o -> (X,omega) -> (
    ZZdfactorization for i to X.period-1 list (-omega^i)*dual((dd^X)_(X.period - i))
    )


--hom
Hom(ZZdFactorization, ZZdFactorization) := ZZdFactorization => (X,Y) -> (
    tensorMF(dual X, Y)
    )

--dHom
Hom(ZZdFactorization, ZZdFactorization, RingElement) := ZZdFactorization => (X,Y,omega) -> (
    dTensor(dual(X,omega),Y,omega)
    )*-

-*OLD
--dual
--note: dual is a method with options
dual(ZZdFactorization) := ZZdFactorization => {} >> o -> X -> (
    ZZdfactorization{-dual((dd^X)_2), dual((dd^X)_1)}
    )

--hom
Hom(ZZdFactorization, ZZdFactorization) := ZZdFactorization => (X,Y) -> (
    tensorMF(dual X, Y)
    )*-



------tensorMF

tensorMF = method()
tensorMF(ZZdFactorization,ZZdFactorization) := (X,Y) -> (
    if not(X.period==2 and Y.period==2) then error "Both inputs should have period 2. Use dTensor method instead.";
    yng := youngest(X,Y);
    if yng.cache#?(tensor,X,Y) then return yng.cache#(tensor,X,Y);
    modules := hashTable for i to 1 list i => (
	directSum for j to 1 list (
	    {j,(i-j)%2} => X_i ** Y_(i-j)
	    )
	);
    A1 := id_(X_0) ** (dd^Y)_1;
    B1 := (dd^X)_1 ** id_(Y_0);	   
    C1 := (dd^X)_2 ** id_(Y_1);
    D1 := -id_(X_1)** (dd^Y)_2;
    diff1 := map(modules#0,modules#1,matrix{{A1, B1},{C1,D1}});
    A2  := id_(X_0) ** (dd^Y)_2;
    B2 := (dd^X)_1 ** id_(Y_1);
    C2 := (dd^X)_2 ** id_(Y_0);
    D2 := -id_(X_1)** (dd^Y)_1;
    diff2 := map(modules#1,modules#0,matrix{{A2, B2},{C2,D2}});
    result := ZZdfactorization{diff1, diff2};
    result.cache.tensor = (X,Y);
    yng.cache#(tensor,X,Y) = result;
    result
)

tensorMF(ZZdFactorizationMap,ZZdFactorizationMap) := (f,g) -> (
    H := new HashTable from {0=>directSum {f_0**g_0,f_1**g_1},1=>directSum {f_0**g_1,f_1**g_0},2 =>directSum {f_0**g_0,f_1**g_1}};
    degf := degree f;
    degg := degree g;
    srcf := source f;
    srcg := source g;
    trgf := target f;
    trgg := target g;
    map(trgf**trgg,srcf**srcg,H,Degree => degf + degg)
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

--i = homological degree
--j = F_j degree
freeMods = method()
freeMods(ZZdFactorization, ZZdFactorization) := (F,G) -> (
    N:= F.period;
for i to N-1 list(
for j to N-1 list(
    {j, (N + i - j)%N} => (F_j)**(G_((N+i-j)%N))
    )))

adjoinRoot = method();
adjoinRoot(ZZ,Ring,Symbol) := (d,Q,t) -> (S1 := Q[t];
    var := (S1_*)_0;
    factored := factor(var^(d)-1);
    cyclo := (factored#(#factored-1))#0;
    S := S1/(cyclo);
    S.rootOfUnity = (S_*)_0;
    S
    )

adjoinRoot(ZZ,Ring,RingElement) := (d,Q,t) -> (adjoinRoot(d,Q,getSymbol "t"))

adjoinRoot(ZZdFactorization,RingElement) := (F,t) -> (
    adjoinRoot(F,getSymbol "t")
    )

adjoinRoot(ZZdFactorization,Symbol) := (F,t) -> (
        d := period F;
	S := adjoinRoot(d,ring F,t);
	P := F**S;
	P.cache.rootOfUnity = (S_*)_0;
	P
	)
    
    adjoinRoot(ZZdFactorizationMap,RingElement) := (F,t) -> (
    adjoinRoot(F,getSymbol "t")
    )

adjoinRoot(ZZdFactorizationMap,Symbol) := (F,t) -> (
        d := period F;
	S := adjoinRoot(d,ring F,t);
	P := F**S;
	P.cache.rootOfUnity = (S_*)_0;
	P
	)
     

--dTensor = method(Options => {RootOfUnity=>false})
-*dTensor = method(Options => {RootOfUnity => true})
dTensor(ZZdFactorization, ZZdFactorization,RingElement) := ZZdFactorization => opts -> (F,G,t) -> (
    --put in check for F.period = G.period
    if (period F)==2 then return tensorMF(F,G);
    if not(opts.RootOfUnity) then return dTensor(F,G,getSymbol "t",RootOfUnity => false);
    Y := youngest(F,G);
    if Y.cache#?(tensor,F,G) then return Y.cache#(tensor,C,D);
    N := F.period;
    modules := hashTable for i to N-1 list i => (
	directSum for j to N-1 list (
	    {j,(j+i)%N} => F_i ** G_(j+i)
	    )
	);
    dF := dd^F;
    dG := dd^G;
    M := for k to N-1 list(
	for i to N-1 list(
	    for j to N-1 list(
		if i == j then (t^i*id_(F_i))**(dG_(k-i+1))		
		else if (j == (i+1)%N) then (dF_(i+1))**(id_(G_((k-i)%N)))
		else 0
	    )));
    result = ZZdfactorization(for i to #M-1 list matrix M_i);
    result.cache.tensor = (F,G);
    Y.cache#(tensor,F,G) = result;
    result
)*-



--dTensor = method(Options => {RootOfUnity=>false})
dTensor = method(Options => {RootOfUnity => true})
dTensor(ZZdFactorization, ZZdFactorization,RingElement) := ZZdFactorization => opts -> (F,G,t) -> (
    --put in check for F.period = G.period
    if (period F)==2 then return tensorMF(F,G);
    if not(opts.RootOfUnity) then return dTensor(F,G,getSymbol "t",RootOfUnity => false);
    Y := youngest(F,G);
    if Y.cache#?(tensor,F,G) then return Y.cache#(tensor,F,G);
    N := F.period;
    modules := hashTable for i to N-1 list i => (
	directSum for j to N-1 list (
	    {j,(i-j)%N} => F_i ** G_(i-j)
	    )
	);
    dF := dd^F;
    dG := dd^G;
    M := for i from 1 to N list(
	map(modules#((i-1)%N),
            modules#(i%N),
            matrix table(
                indices modules#((i-1)%N),
                indices modules#(i%N),
                (j,k) -> (
                    tar := component(modules#((i-1)%N), j);
                    src := component(modules#(i%N), k);
                    map(tar, src, 
                        if ({(k#0-j#0)%N,(k#1-j#1)%N} == {0,1}) then (t^(k#0)*id_(F_(k#0)))**(dG_(k#1))
                        else if ({(k#0-j#0)%N,(k#1-j#1)%N} == {1,0})  then (dF_(k#0))**(id_(G_(k#1)))
                        else 0)
                    ))));
    result := ZZdfactorization(M);
    result.cache.tensor = (F,G);
    Y.cache#(tensor,F,G) = result;
    result
)



dTensor(ZZdFactorization,ZZdFactorization,Symbol) := ZZdFactorization => opts -> (F,G,t) -> (S := adjoinRoot(period F,ring F,t);
    dTensor(F**S,G**S,(S_*)_0)
    )

dTensor(List,RingElement) := ZZdFactorization => opts -> (L,t) -> (
    if #L==2 then return dTensor(L_0,L_1,t,opts);
    if not(opts.RootOfUnity) then return dTensor(L,getSymbol "t");
    dTensor(dTensor(L_{0..#L-2},t),L_(#L-1),t)
    )

dTensor(List,Symbol) := ZZdFactorization => opts -> (L,t) -> (
    if #L==2 then return dTensor(L_0,L_1,t,opts);
    S := adjoinRoot(period L_0,ring L_0,t);
    Ln := apply(L,i->i**S);
    dTensor(dTensor(Ln_{0..#Ln-2},(S_*)_0),Ln_(#Ln-1),(S_*)_0)
    )

dTensor(ZZdFactorizationMap,ZZdFactorizationMap,RingElement) := ZZdFactorizationMap => opts -> (f,g,t) -> (
    if not(opts.RootOfUnity) then return dTensor(f,g,getSymbol "t");
    degf := degree f;
    degg := degree g;
    srcf := source f;
    srcg := source g;
    trgf := target f;
    trgg := target g;
    d := period srcf;
    L = for i to d list (
	     for j to d-1 list (j,(j+i)%d)
	     );
    Ln = new HashTable from toList apply(0..d,i->i=>directSum apply(L_i,j->(f_(j_0)**g_(j_1))));
    map(trgf**trgg,srcf**srcg,Ln,Degree => degf + degg)
    )

dTensor(ZZdFactorizationMap,ZZdFactorizationMap,Symbol) := ZZdFactorizationMap => opts -> (f,g,t) -> (
    S := adjoinRoot(period source f,ring f,t);
    dTensor(f**S,g**S,(S_*)_0)
    )
	    
    
tensor(ZZdFactorization,ZZdFactorization) := ZZdFactorization => {} >> opts -> (F,G) -> (
    if not(F.period == G.period) then error "Expected factorizations with the same period";
    if F.period==2 then return tensorMF(F,G);
    if (ring F).?rootOfUnity then return dTensor(F,G,(ring F).rootOfUnity)
    else error "Must adjoin dth root of unity when input has period d > 2";
    )

tensor(ZZdFactorization,ZZdFactorization,RingElement) := ZZdFactorization => {Dispatch => {ZZdFactorization,ZZdFactorization,RingElement}} >> opts -> (F,G,t) -> (print("here we");
    if not(F.period==G.period) then error "Expected factorizations with the same period";
    if F.period==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dTensor(F,G,t,RootOfUnity=>false)
    )

tensor(ZZdFactorization,ZZdFactorization,Symbol) := ZZdFactorization => {} >> opts -> (F,G,t) -> (
    if not(F.period==G.period) then error "Expected factorizations with the same period";
    if F.period==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dTensor(F,G,t)
    )

tensor(ZZdFactorizationMap,ZZdFactorizationMap) := ZZdFactorizationMap => {} >> opts -> (f,g) -> (
    F := source f;
    G := source g;
    if not(F.period == G.period) then error "Expected factorizations with the same period";
    if F.period==2 then return tensorMF(f,g);
    if (ring F).?rootOfUnity then return dTensor(f,g,(ring F).rootOfUnity)
    else error "Must adjoin dth root of unity when input has period d > 2";
    )

tensor(ZZdFactorizationMap,ZZdFactorizationMap,RingElement) := ZZdFactorizationMap => {} >> opts -> (f,g,t) -> (
    F := source f;
    G := source g;
    if not(F.period==G.period) then error "Expected factorizations with the same period";
    if F.period==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dTensor(f,g,t,RootOfUnity=>false)
    )

tensor(ZZdFactorizationMap,ZZdFactorizationMap,Symbol) := ZZdFactorizationMap => {} >> opts -> (f,g,t) -> (
    F := source f;
    G := source g;
    if not(F.period==G.period) then error "Expected factorizations with the same period";
    if F.period==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dTensor(f,g,t)
    )


ZZdFactorization ** ZZdFactorization := ZZdFactorization => (F,G) -> (tensor(F,G))
Complex ** ZZdFactorization := ZZdFactorization => (C,F) -> (tensor(Fold(C**(ring F),period F),F))
ZZdFactorization ** Complex := ZZdFactorization => (F,C) -> (tensor(F,Fold(C**(ring F),period F)))
ZZdFactorization ** ZZdFactorizationMap := ZZdFactorizationMap => (F,f) -> (tensor(id_F,f))
ZZdFactorizationMap ** ZZdFactorization := ZZdFactorizationMap => (f,F) -> (tensor(f,id_F))
ComplexMap ** ZZdFactorizationMap := ZZdFactorizationMap => (f,g) -> (tensor(Fold(f**(ring g),period source g),g))
ZZdFactorizationMap ** ComplexMap := ZZdFactorizationMap => (f,g) -> (tensor(f,Fold(g**(ring f),period source f)))
ZZdFactorizationMap ** ZZdFactorizationMap := ZZdFactorizationMap => (f,g) -> (tensor(f,g))
 





dual(ZZdFactorization) := ZZdFactorization => {} >> opts -> F -> (
    if F.period == 2 then return ZZdfactorization {-dual F.dd_2,dual F.dd_1};
    if (ring F).?rootOfUnity then return dDual(F,(ring F).rootOfUnity)
    else error "Must adjoin dth root of unity when input has period d > 2";
    )

dual(ZZdFactorization,RingElement) := ZZdFactorization => {} >> opts -> (F,t) -> (
    if F.period ==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dDual(F,t,RootOfUnity=>false)
    )

dual(ZZdFactorization,Symbol) := ZZdFactorization => {} >> opts -> (F,t) -> (
    if F.period ==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dDual(F,t)
    )

dual(ZZdFactorizationMap) := ZZdFactorizationMap => {} >> opts -> f -> (
    deg := degree f;
    F = source f;
    G = target f;
    if F.period==2 then return map(dual F,dual G,new HashTable from {1=>(dual f_0),2=>(dual f_1)},Degree => -degree f);
    if (ring F).?rootOfUnity then return dDual(f,(ring F).rootOfUnity)
    else error "Must adjoin dth root of unity when input has period d > 2";
    )

dual(ZZdFactorizationMap,RingElement) := ZZdFactorizationMap => {} >> opts -> (f,t) -> (
    if F.period ==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dDual(f,t,RootOfUnity=>false)
    )

dual(ZZdFactorizationMap,Symbol) := ZZdFactorizationMap => {} >> opts -> (f,t) -> (
    if F.period ==2 then error "No need to specify root of unity for ZZ/2-graded factorization";
    dDual(f,t)
    )

dDual = method(Options => {RootOfUnity => true});
dDual(ZZdFactorization,RingElement) := ZZdFactorization => opts -> (F,t) -> (
    if not(opts.RootOfUnity) then return dDual(F,getSymbol "t");
    diffs := reverse values (F.dd.map);
    ZZdfactorization toList apply(0..#diffs-1,i->-t^i*dual(diffs#i))
    )

dDual(ZZdFactorization,Symbol) := ZZdFactorization => opts -> (F,t) -> (
    S := adjoinRoot(period F,ring F,t);
    dDual(F**S,(S_*)_0)
    )

dDual(ZZdFactorizationMap,RingElement) := ZZdFactorization => opts -> (f,t) -> (
    if not(opts.RootOfUnity) then return dDual(t,getSymbol "t");
    dualMaps := reverse values (f.map);
    d := #dualMaps;
    srcf := source f;
    trgf := target f;
    degf := degree f;
    map(dual srcf, dual trgf, new HashTable from toList apply(1..d,i->i=>dual(dualMaps#(i%d))),Degree => -degf)
    )

dShift =  method(Options => {RootOfUnity => true});
dShift(ZZ,ZZdFactorization,RingElement) := ZZdFactorization => opts -> (s,F,t) -> (
    if not(opts.RootOfUnity) then return dShift(s,F,getSymbol "t");
    d := period F;
    diffs := values (F.dd.map);
    ZZdfactorization toList apply(0..#diffs-1,i->t^s*diffs#((i+s)%d))
    )

dShift(ZZdFactorization,Symbol) := ZZdFactorization => opts -> (s,F,t) -> (
    S := adjoinRoot(period F,ring F,t);
    dDual(s,F**S,(S_*)_0)
    )

ZZdFactorization Array := (C, L) -> (
    if #L != 1 or not instance(L#0,ZZ) then error "expected an integer shift";
    if period C == 2 then return ZZdfactorization {(-1)^(L#0)*C.dd_(L#0+1),(-1)^(L#0)*C.dd_(L#0)};
    if (ring C).?rootOfUnity then return dShift(L#0,C,(ring C).rootOfUnity)
    else error "Must adjoin dth root of unity when input has period d > 2";
    )


    

--this gives a trivial d-fold factorization of a monomial
trivialFactorization = method()
trivialFactorization(RingElement) := f -> (if not(#(terms f)==1) then error "Expected ring element to be monomial";
    Lk := toList factor(f);
    theDiffs = (flatten(Lk/(i->toList(i#1:i#0))))/(j->matrix{{j}});
    if first degree f < #theDiffs then theDiffs = {(last theDiffs)*first(theDiffs)}|(theDiffs_{1..#theDiffs-2});
    ZZdfactorization theDiffs
    )

linearFactorization = method(Options => {RootOfUnity => true});
linearFactorization(RingElement) := ZZdFactorization => opts -> (f) -> (
    if not isHomogeneous f then error "Expected homogeneous element";
    if not((ring f).?rootOfUnity) and not first degree(f)==2 then error "Must adjoint dth root of unity if degree > 2";
    L := (terms f)/trivialFactorization;
    if first degree f == 2 then return tensor(L);
    if not(opts.RootOfUnity) then return dTensor(L,getSymbol "t",opts);
    dTensor(L,t)
    )

linearFactorization(RingElement,RingElement) := ZZdFactorization => opts -> (f,t) -> (
    linearFactorization(f,getSymbol "t")
    )

linearFactorization(RingElement,Symbol) := ZZdFactorization => opts -> (f,t) -> (
    S := adjoinRoot(first degree f,ring f,t);
    linearFactorization(sub(f,S))
    )

randomFactorization = method();
randomFactorization(ZZ,Ring) := (d,Q) -> (
    f := random(d,Q);
    linearFactorization(f)
    )

randomFactorization(ZZ,Ring,RingElement) := (d,Q,t) -> (
    f := random(d,Q);
    linearFactorization(f,t)
    )

randomFactorization(ZZ,Ring,Symbol) := (d,Q,t) -> (
    f := random(d,Q);
    linearFactorization(f,t)
    )

end--


restart
needs "MF_functions.m2"

S = QQ[x,y,z,u,v,w]
P = ZZdfactorization {x,y,z,u,v,w}
F = ZZdfactorization{matrix{{x,0},{0,x}}, matrix{{y,0},{0,y}}, matrix{{z,0},{0,z}}}
G = ZZdfactorization{matrix{{v,x},{0,v}}, matrix{{v,x},{x,v}}, matrix{{w,0},{0,w}}}

isdFactorization F --true
isdFactorization G --false

--put in check for tensorMF
Q = adjoinRoot(F.period, S)
T = dTensor(F,G,t)
T.period
dd^T
(dd^T)_0*(dd^T)_1*(dd^T)_2
isdFactorization T --false; do we want it false?


---more test
R = QQ[x,y,z,w]
K = koszulMF({x^2*y^2,z^2*w^2, x*y, y^10, x*y, z^100})
KMF = koszulMF({x,x,y,y,z,z})
isdFactorization KMF

----test
R = QQ[a,b,c,d]
K = koszulMFf({a,b,c,d^2}, a^5*d + d^100 + a*b*c)
isMatrixFactorization K
dd^K
I = ideal {a,b,c,d^2}
f = a^5*d + d^100 + a*b*c
koszulMFf(I, f)

M = f //gens(I)
N = transpose (gens I) | M
E = entries N
F = select(E, i -> not(i#1 == 0))
flatten F

K = koszulMFf({a,b,c,d}, a^3 + b^3 +c^3 +d^3)
dd^K

--test with larger tensor
restart
needs "MF_functions.m2"
Q=QQ[x_1..x_3];
F = ZZdfactorization {matrix{{x_1}},matrix{{x_2}},matrix{{x_3}}}; -- defining two trivial factorizations
G = ZZdfactorization {matrix{{x_2}},matrix{{x_2}},matrix{{x_3}}};
S = adjoinRoot(3,Q)
peek S
T = dTensor(F**S,G**S,t)
T.dd
F.dd
(dd^T_1)*(dd^T_2)*(dd^T_3)

R = QQ[x,y,z,u,v,w]
Q = adjoinRoot(4,R)
koszulMF{{x,y,z},{u,v,w},{x^10*y, z,w},{u^2,v,w}}

koszulMF = method()
koszulMF(List, RingElement) := (L, omega) -> (
    n := #(L#0);
    d := #L;
    Z := for i to n-1 list ZZdfactorization for j to d-1 list (L#j)#i;
    T := Z#0;
    for i from 1 to n-1 do T = dTensor(T,Z#i, omega);
    T
)

K = koszulMF({{x,y,z},{u,v,w},{x^10*y, z,w},{u^2,v,w}}, t)

dd^K

isdFactorization K

koszulMFf(List, RingElement, RingElement) := (L,f, omega) -> (
    I := ideal L;
    M := f //gens(I^2);
    N := transpose (gens I) | M;     
    E  := entries N;
    F := select(E, i -> not(i#1 == 0)); 
    A := flatten F;
    koszulMF(A)
)

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
	 
S = QQ[x,y,z,w]
Q = adjoinRoot(5, S)
K = koszulMFf({x,y,z,w}, x^5+y^5+w^5, 4, t)
dd^K
isdFactorization K

