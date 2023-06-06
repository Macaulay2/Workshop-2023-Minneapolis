--restart
load "ZZdFactorizations.m2" 

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

KoszulMF = method()
KoszulMF(List) := L -> (

X := ZZdfactorization{ matrix{{L#0}}, matrix{{L#1}}};
for i from 1 to (#L)//2-1 do X = tensorMF(X, ZZdfactorization{ matrix{{L#(2*i)}}, matrix{{L#(2*i+1)}}});

X
)


------KoszulMF from ideal and function
KoszulMFf = method()
KoszulMFf(List, RingElement) := (L,f) -> (
    I := ideal L;
    M := f //gens(I);
    N := transpose (gens I) | M;     
    E  := entries N;
    F := select(E, i -> not(i#1 == 0)); 
    A := flatten F;
    KoszulMF(A)
)

KoszulMFf(Ideal, RingElement) := (I,f) -> (
    KoszulMFf(flatten entries gens I, f)
)

eulerMF = f -> (KoszulMFf(ideal jacobian f , f))


end--

needs "tensorMF.m2"

---test
R = QQ[x,y,z,w]
K = KoszulMF({x^2*y^2,z^2*w^2, x*y, y^10, x*y, z^100})

KMF = KoszulMF({x,x,y,y,z,z})

(dd^KMF)_1*(dd^KMF)_2

----test
R= QQ[a,b,c,d]
K = KoszulMFf({a,b,c,d^2}, a^5*d + d^100 + a*b*c)
isMatrixFactorization K
dd^K
I = ideal {a,b,c,d^2}
f = a^5*d + d^100 + a*b*c
KoszulMFf(I, f)

M = f //gens(I)
N = transpose (gens I) | M
E = entries N
F = select(E, i -> not(i#1 == 0))
flatten F

K = KoszulMFf({a,b,c,d}, a^3 + b^3 +c^3 +d^3)
dd^K

K= KoszulMFf({a,b,c}, a^2-b^2+c^5)

K= KoszulMFf({a+b,c}, a^2-b^2+c^5)

K= KoszulMFf({a+b,c^3}, a^2-b^2+c^5)
