------tensor product

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
