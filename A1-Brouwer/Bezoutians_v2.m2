--Nikita, Gabriel
-- Given an endomorphism of affine space f=(f1,...,fn) given as a list of polynomials, return the Bezoutian corresponding to the endomorphism
bezoutian = method()

bezoutian (List) := (Matrix) => (Endo) -> (
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
    mapxtoX:= (map(RX,S,toList(X_1..X_n))); -- defines the map f_i(x_1, ..., x_n) to f_i(X_1, ..., X_n)
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n))); -- defines the map f_i(x_1, ..., x_n) to f_i(Y_1, ..., Y_n)
    standBasisX := basis (RX/(ideal (leadTerm (mapxtoX ideal Endo)))); -- Compute the standard basis of kk[X_1, ..., X_n]/(f_1, ..., f_n)
    standBasisY := basis (RY/(ideal (leadTerm (mapxtoY ideal Endo)))); -- Compute the standard basis of kk[Y_1, ..., Y_n]/(f_1, ..., f_n)
    id1 := (ideal apply(toList(0..n-1), i-> mapxtoX(Endo_i))); -- defines an ideal (f_1(X), ..., f_n(X))
    id2 := (ideal apply(toList(0..n-1), i-> mapxtoY(Endo_i))); -- defines an ideal (f_1(Y), ..., f_n(Y))
    promotedEndo :=sub(id1,R)+sub(id2,R); -- takes the sum of the ideals in the ring kk[X_1..Y_n]
    Rquot:= R/promotedEndo; -- defines the quotient ring kk[X_1..Y_n]/(f(X),f(Y)) which is Q(f)\otimes_{k}Q(f) in the paper
    sBXProm :=sub(standBasisX,Rquot); -- moves the standard bases to the quotient ring
    sBYProm :=sub(standBasisY,Rquot); -- moves the standard bases to the quotient ring 
    bezDetProm := promote(bezDet, Rquot); -- moves the Bezoutian polynomial to the quotient ring
    ------------------------
    phi0 := map(kk,Rquot,(toList ((2*n):0))); -- define a ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring (RMK is this even needed?)
    --will return matrix B
    m:= numColumns(sBXProm);
    B:= mutableMatrix id_(kk^m);
    --print sBXProm#1;
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            B_(i,j)=phi0(coefficient((sBXProm_(0,i)**sBYProm_(0,j))_(0,0), bezDetProm));
        );
    );
    return matrix(B);
);

T = QQ[x_1,x_2];
e = {x_1^2+x_2,x_2^3+x_1};
bez_e = bezoutian(e);
bez_e = transpose bez_e;
print bez_e


T1 = QQ[x_1];
f = {x_1^4 + x_1^3 - x_1^2 - x_1};
bez_f = bezoutian(f);
bez_f = transpose bez_f;
print bez_f



T = QQ[x_1,x_2];
e = {x_1^2+x_2,x_2^3+x_1};
print bezoutian(e);
-- Given an endomorphism of an affine space f=(f1,...,fn) given given as a sequence, and a closed point p of the target, return the local Bezoutian
--localBezoutian = method()
