
newPackage(
    "examplePackage",
    DebuggingMode=>true
    )

export{
    "myFunction"
    }


myFunction = method()
myFunction(List) := ZZ => ( L )->(
    n := #L;
    kk := coefficientRing(ring(L#0));
    X:= symbol X;
    R := kk[X_1..X_n];
    return #gens(R);
    

);
S = QQ[x,y];
L1 = {x^2,x*y};
myFunction(L1)
