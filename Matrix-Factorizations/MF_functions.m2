needs "ZZdFactorizations.m2" 

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


--dual
--note: dual is a method with options
dual(ZZdFactorization) := ZZdFactorization => {} >> o -> X -> (
    ZZdfactorization{-dual((dd^X)_2), dual((dd^X)_1)}
    )


--hom
Hom(ZZdFactorization, ZZdFactorization) := ZZdFactorization => (X,Y) -> (
    tensorMF(dual X, Y)
    )

------tensorMF

tensorMF = method()
tensorMF(ZZdFactorization,ZZdFactorization) := (X,Y) -> (
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


-----KoszulMF
--koszul matrix factorization
--in: list of even number of elements in ring
--out: matrix factorization of sum products of pairs (L#i)*(L#(i+1))


koszulMF = method()
koszulMF(List) := L -> (
    X := ZZdfactorization{ matrix{{L#0}}, matrix{{L#1}}};
    for i from 1 to (#L)//2-1 do X = tensorMF(X, ZZdfactorization{ matrix{{L#(2*i)}}, matrix{{L#(2*i+1)}}});
    X
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

end--


restart
needs "MF_functions.m2"

--building N-fold tensor

S = QQ[x,y,z, u, v, w]
F = ZZdfactorization{matrix{{x}}, matrix{{y}}, matrix{{z}}}
G = ZZdfactorization{matrix{{u}}, matrix{{v}}, matrix{{w}}}
(dd^F)_(-50)
dd^F
--add check that F.period = G.period
N = F.period

--i = homological degree
--j = F_j degree
freeMods = method()
freeMods(ZZdFactorization, ZZdFactorization) := (F,G) -> (
    N:= F.period;
for i to N-1 list(
for j to N-1 list(
    {j, (N + i - j)%N} => (F_j)**(G_((N+i-j)%N))
    )))

freeMods(F, G)
netList oo

dTensor = method()
dTensor(ZZdFactorization, ZZdFactorization) := (F,G) -> (
    N := F.period;
    dF := dd^F;
    dG := dd^G;
    for i to N-1 list(
	for j to N-1 list(
	    {j, (N + i - j)%N} =>
	    if (i==j) (id_(F_i)**(dG)_(j+1))
	    else if (j==((i+1)%N))
	    else 0 
    
k=0
N = F.period
dF = dd^F
dG = dd^G
for i to N-1 list(
	for j to N-1 list(
	    {j, (i)%N} =>(
	    if i==j then id_(F_i)**(dG)_(k-i+1)
	    --else if then 0 --(j==((i+1)%N))
	    else 0)

)
)
netList oo


M = for k to N-1 list(
    for i to N-1 list(
	for j to N-1 list(
	    if i == j then 
		(id_(F_i))**(dG_((k-i+1)%N))		
	    --else if (j == (i+1)%N) then(
		--0
		--)
	    else 0
	    ))
    )
netList oo
matrix M_0
M_0
id_(F_2)**(dG_((0-2+1)%3))

--if i < 0 then "neg" 
 --              else if i == 0 then "zer" 
   --            else "pos")
--tests


---test
R = QQ[x,y,z,w]
K = koszulMF({x^2*y^2,z^2*w^2, x*y, y^10, x*y, z^100})

KMF = koszulMF({x,x,y,y,z,z})

(dd^KMF)_1*(dd^KMF)_2

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
