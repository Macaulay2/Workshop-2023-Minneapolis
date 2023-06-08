
-- Given an endomorphism of affine space f=(f1,...,fn) given as a list of polynomials, return the Bezoutian corresponding to the endomorphism
--for notation: input domain is kk[x_1..x_n] while Bezoutian computed in kk[X_1..Y_n]

load "localAlgebraBasis.m2"

localA1Degree = method()

localA1Degree (List, Ideal) := (Matrix) => (Endo,p) -> (
    -- Endo is the list {f_1, f_2, ..., f_n} of polynomials kk^n -> kk^n
     n := #Endo; -- n is the no. of polynomials defining the morphism    
    kk := coefficientRing(ring(Endo#0)); -- kk is the field, where f is an endomorphsism between affine n-space over kk. 
     S:=ring(Endo#0); -- S is the ring in which the defining polynomials of f live
    
    -- first check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then (print "Error: morphism does not have isolated zeroes"; return Endo;);
    
    -- create internal rings/matrices
   
    R := kk[X_1..Y_n]; -- inititalize this new ring for the computation of the Bezoutian.
    D := mutableMatrix id_((frac R)^n); -- Create a (n x n) matrix D which will be populated by \Delta_{ij} in the paper
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(
	    -- iterate through the entries of the matrix D and populate it with the following information ...
	    targetList1 := toList (Y_1..Y_j|X_(j+1)..X_n); -- create the list {x_1 => Y_1, ..., x_j-1 => Y_j-1, x_j => X_j, ..., x_n => X_n}
	    targetList2 := toList (Y_1..Y_(j+1)| X_(j+2)..X_n); -- create the list {x_1 => Y_1, ..., x_j => Y_j, x_j+1 => X_j+1, ..., x_n => X_n}
        numeratorD = ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i)); -- map f_i(x) to f_i(Y_1, ..., Y_j-1, X_j, ..., X_n) resp. 
	    D_(i,j)= numeratorD/(X_(j+1)-Y_(j+1)); -- this is \Delta_{i,j} from the paper
	); 
    );

    -- bezDet is the determinant of the (n x n) matrix D
    bezDet:= numerator (det(D)); -- typecast to kk[X_1, ..., Y_n]

    -- define formal variabls X_i, Y_i that replace x_i
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

    -- Create the ideal J.  It has prop that  k[x_1 .. x_n]_p/I_p is isom to k[x_1 .. x_n]/J
    J:=(ideal Endo):saturate(ideal Endo,p);
    
    localIdeal :=sub(mapxtoX(J),R)+sub( mapxtoY(J),R); -- create copies of J with the formal variables
    Rquot:=R/localIdeal;  -- this is the tensor product  Q_p(f) x  Q_p(f), where X's are the variables in first term, Y's are variables in second part;
    
    sBXProm :=apply(toList(0..#standBasisX-1),i-> sub(standBasisX_i,Rquot)); -- moves the standard bases to the quotient ring
    sBYProm :=apply(toList(0..#standBasisY-1),i-> sub(standBasisY_i,Rquot)); -- moves the standard bases to the quotient ring 
   
    -- Now reduce the bezDet determinant subject ot the local ideal in both the X's, Y's.
    bezDetRed := bezDet % localIdeal;
   
    
    ------------------------
    
    phi0 := map(kk,Rquot,(toList ((2*n):0))); -- ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring 
    
    -- m is the dimension of the basis for the local ring
    m:= #sBXProm;
    
    -- Now create Bezoutian matrix B for the quadratic form by reading off the local coefficients. 
    -- B is a (m x m) matrix.  Coefficent B_(i,j) is the coefficient of the (ith basis vector x jth basis vector) in tensor product.
    -- phi0 maps the coefficient to kk
    
    B:= mutableMatrix id_(kk^m);
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            B_(i,j)=phi0(coefficient((sBXProm_i**sBYProm_j)_(0,0), bezDetRed));
        );
    );
    return matrix(B);
);


