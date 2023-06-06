restart
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
for i from 1 to (#L)//2-1 do X = tensorMF(X, ZZdfactorization{ L#(2*i)), L#(2*i+1)});

X
)

---test
R = QQ[x,y,z,w]
K = KoszulMF({x^2*y^2,z^2*w^2, x*y, y^10, x*y, z^100})

KMF = KoszulMF({x,x,y,y,z,z})

(dd^KMF)_1*(dd^KMF)_2
