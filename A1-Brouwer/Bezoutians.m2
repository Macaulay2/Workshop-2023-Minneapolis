-- Given an endomorphism of affine space f=(f1,...,fn) given as a list of polynomials, return the Bezoutian corresponding to the endomorphism
bezoutian = method()

bezoutian (List) := (Matrix) => (Endo) -> (
    --use ring ideal(Endo);
    if dim ideal(Endo) > 0  then (print "Error: morphism does not have isolated zeroes"; return Endo;);
    n := #Endo;
    kk := coefficientRing(ring(Endo#0));
    S:=ring(Endo#0);
    R := kk[X_1..Y_n]; -- Here need to initialize R as kk[X_1..Y_n,x_1..x_n]???? Can substitute and change ring simultaneously? 
    D := mutableMatrix id_((frac R)^n);
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(
	    -- substitute the variables and divide by the X_i - Y_j
	    -- Need to find a way to combine two lists in the following way: eg. {1,2,3,4,5} and {a,b,c,d,e} go to {1=>a, 2=>b, etc.}
	    targetList1 := toList (Y_1..Y_j|X_(j+1)..X_n);
	    targetList2 := toList (Y_1..Y_(j+1)| X_(j+2)..X_n);
        numeratorD = ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i));
	    D_(i,j)= numeratorD/(X_(j+1)-Y_(j+1));
	); 
    );
    bezDet:= numerator (det(D));--take numerator to be polynomial (not rational fxn)
    RX:=kk[X_1..X_n];
    RY:=kk[Y_1..Y_n];
    mapxtoX:= (map(RX,S,toList(X_1..X_n)));
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n)));
    standBasisX := basis (RX/(ideal (leadTerm (mapxtoX ideal Endo))));
    standBasisY := basis (RY/(ideal (leadTerm (mapxtoY ideal Endo))));
    id1 := (ideal apply(toList(0..n-1), i-> mapxtoX(Endo_i)));
    id2 := (ideal apply(toList(0..n-1), i-> mapxtoY(Endo_i)));
    promotedEndo :=sub(id1,R)+sub(id2,R);
    Rquot:= R/promotedEndo;
    sBXProm :=sub(standBasisX,Rquot);
    sBYProm :=sub(standBasisY,Rquot);
    --sBRquot := sBXProm**sBYProm;
    bezDetProm := promote(bezDet, Rquot);
    ------------------------
    phi0 := map(kk,Rquot,(toList ((2*n):0)));
    coeffsTuple := coefficients(bezDetProm);
    presMons :=(coeffsTuple)_0;  
    coeffs := (coeffsTuple)_1;
    print class phi0((coeffs)_(2,0));
    --will return matrix B
    m:=#sBXProm;
    B:=mutableMatrix id_(kk^m);
    coefIndex:=0;
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            --sBXProm_
            B_(i,j)=phi0((coeffs)_(coefIndex,0));
        );
    );
    return 1;
);

T = QQ[x_1,x_2];
e = {x_1^2+x_2,x_2^3+x_1};
print bezoutian(e);
-- Given an endomorphism of an affine space f=(f1,...,fn) given given as a sequence, and a closed point p of the target, return the local Bezoutian
--localBezoutian = method()
