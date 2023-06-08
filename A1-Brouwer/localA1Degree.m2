-- Given an endomorphism of affine space f=(f1,...,fn) given as a list of polynomials, return the Bezoutian corresponding to the endomorphism
--for notation: input domain is kk[x_1..x_n] while Bezoutian computed in kk[X_1..Y_n]
load "localAlgebraBasis.m2"

localA1Degree = method()

localA1Degree (List, Ideal) := (Matrix) => (Endo,p) -> (
    -- first check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then (print "Error: morphism does not have isolated zeroes"; return Endo;);
    n := #Endo; -- n is the no. of polynomials defining the morphism
    kk := coefficientRing(ring(Endo#0)); -- kk is the field, where f is an endomorphsism between affine n-space over kk. 
    S:=ring(Endo#0); -- S is the ring in which the defining polynomials of f live
    R := kk[X_1..Y_n]; -- inititalize this new ring for the computation of the Bezoutian.
    D := mutableMatrix id_((frac R)^n); -- Create a matrix which will be populated by \Delta_{ij} in the paper
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(
	    -- iterate through the entries of the matrix D and populate it with the following information ...
	    targetList1 := toList (Y_1..Y_j|X_(j+1)..X_n); -- create the list {x_1 => Y_1, ..., x_j-1 => Y_j-1, x_j => X_j, ..., x_n => X_n}
	    targetList2 := toList (Y_1..Y_(j+1)| X_(j+2)..X_n); -- create the list {x_1 => Y_1, ..., x_j => Y_j, x_j+1 => X_j+1, ..., x_n => X_n}
        numeratorD = ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i)); -- map f_i(x) to f_i(Y_1, ..., Y_j-1, X_j, ..., X_n) resp. 
	    D_(i,j)= numeratorD/(X_(j+1)-Y_(j+1)); -- this is \Delta_{i,j} from the paper
	); 
    );
    bezDet:= numerator (det(D)); -- typecast to kk[X_1, ..., Y_n]
    RX:=kk[X_1..X_n]; 
    RY:=kk[Y_1..Y_n];
    mapxtoX:= (map(RX,S,toList(X_1..X_n))); -- map f_i(x_1, ..., x_n) to f_i(X_1, ..., X_n)
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n))); -- map f_i(x_1, ..., x_n) to f_i(Y_1, ..., Y_n)

    --find standard basis and define local quotient ring
    list1 := (apply(toList(0..n-1), i-> mapxtoX(Endo_i))); -- list (f_1(X), ..., f_n(X))
    list2 := (apply(toList(0..n-1), i-> mapxtoY(Endo_i))); -- list (f_1(Y), ..., f_n(Y))
    --apply localAlgebraBasis to compute standard basis for localization
    standBasisX :=localAlgebraBasis(list1,mapxtoX p);
    standBasisY := localAlgebraBasis(list2,mapxtoY p); 
    promotedEndo :=sub(ideal list1,R)+sub(ideal list2,R); -- takes the sum of the ideals in the ring kk[X_1..Y_n]
    Rquot:= R/promotedEndo; -- localized quotient ring
    
    sBXProm :=apply(toList(0..#standBasisX-1),i-> sub(standBasisX_i,Rquot)); -- moves the standard bases to the quotient ring
    sBYProm :=apply(toList(0..#standBasisY-1),i-> sub(standBasisY_i,Rquot)); -- moves the standard bases to the quotient ring 
    bezDetProm := promote(bezDet, Rquot); -- moves the Bezoutian polynomial to the quotient ring
    ------------------------
    phi0 := map(kk,Rquot,(toList ((2*n):0))); -- ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring (RMK is this even needed?)
    --will return matrix B
    m:= #sBXProm;
    B:= mutableMatrix id_(kk^m);
    --print sBXProm#1;
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            B_(i,j)=phi0(coefficient((sBXProm_i**sBYProm_j)_(0,0), bezDetProm));
        );
    );
    return matrix(B);
);

T = QQ[z];
e = {z^4+z^3-z^2-z};
p=ideal {z+1/1};
bez_e = localA1Degree(e,p);
bez_e = transpose bez_e;
print bez_e