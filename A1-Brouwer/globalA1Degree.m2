--Nikita, Gabriel
-- Given an endomorphism of affine space f=(f1,...,fn) given as a list of polynomials, return the Bezoutian corresponding to the endomorphism
-- For notation: input domain is kk[x_1..x_n] while Bezoutian computed in kk[X_1..Y_n]


globalA1Degree = method()

globalA1Degree (List) := (Matrix) => (Endo) -> (
    -- Endo is the list {f_i} of the polynomials f_i
    -- Check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then (print "Error: morphism does not have isolated zeroes"; return Endo;);
    
    -- Number of polynomials defining the morphism
    n := #Endo;
    
    -- The morphism f is an endomorphsism between affine n-space over field kk 
    kk := coefficientRing(ring(Endo#0));
    
    -- S is the ring in which the defining polynomials of f live
    S:=ring(Endo#0);
    
    -- Inititalize this new ring for the computation of the Bezoutian
    R := kk[X_1..Y_n];
    
    -- Create a matrix which will be populated by \Delta_{ij} in the paper
    D := mutableMatrix id_((frac R)^n);
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(
	    -- Iterate through the entries of the matrix D and populate it with the following
	    
	    -- Create the list {x_1 => Y_1, ..., x_j-1 => Y_j-1, x_j => X_j, ..., x_n => X_n}
	    targetList1 := toList (Y_1..Y_j|X_(j+1)..X_n);
	    
	    -- Create the list {x_1 => Y_1, ..., x_j => Y_j, x_j+1 => X_j+1, ..., x_n => X_n}
	    targetList2 := toList (Y_1..Y_(j+1)| X_(j+2)..X_n);
	    
	-- Map f_i(x) to f_i(Y_1, ..., Y_j-1, X_j, ..., X_n) resp. 
        numeratorD = ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i));
	
	    -- This is \Delta_{i,j} from the paper
	    D_(i,j)= numeratorD/(X_(j+1)-Y_(j+1)); 
	); 
    );

    -- Typecast to kk[X_1, ..., Y_n]
    bezDet:= numerator (det(D));
    RX:=kk[X_1..X_n]; 
    RY:=kk[Y_1..Y_n];
    
    -- Defines the map f_i(x_1, ..., x_n) to f_i(X_1, ..., X_n)
    mapxtoX:= (map(RX,S,toList(X_1..X_n)));
    
    -- Defines the map f_i(x_1, ..., x_n) to f_i(Y_1, ..., Y_n)
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n)));
    
    -- Compute the standard basis of kk[X_1, ..., X_n]/(f_1, ..., f_n)
    standBasisX := basis (RX/(ideal (leadTerm (mapxtoX ideal Endo))));
    
    -- Compute the standard basis of kk[Y_1, ..., Y_n]/(f_1, ..., f_n)
    standBasisY := basis (RY/(ideal (leadTerm (mapxtoY ideal Endo))));
    
    -- Defines an ideal (f_1(X), ..., f_n(X))
    id1 := (ideal apply(toList(0..n-1), i-> mapxtoX(Endo_i)));
    
    -- Defines an ideal (f_1(Y), ..., f_n(Y)) 
    id2 := (ideal apply(toList(0..n-1), i-> mapxtoY(Endo_i)));
    
    -- Takes the sum of the ideals in the ring kk[X_1..Y_n]
    promotedEndo :=sub(id1,R)+sub(id2,R); 
    
    -- Defines the quotient ring kk[X_1..Y_n]/(f(X),f(Y)) which is Q(f)\otimes_{k}Q(f) in the paper
    Rquot:= R/promotedEndo;
    
    -- Moves the standard bases to the quotient ring
    sBXProm :=sub(standBasisX,Rquot); 
    
    -- Moves the standard bases to the quotient ring
    sBYProm :=sub(standBasisY,Rquot);
    
    --Reduces bezDet mod promotedEndo 
    bezDetRed := bezDet % promotedEndo;
    ------------------------
    
    -- Define a ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring (RMK is this even needed?)
    phi0 := map(kk,Rquot,(toList ((2*n):0)));
    -- Returns matrix B
    m:= numColumns(sBXProm);
    B:= mutableMatrix id_(kk^m);
    
    -- Print sBXProm#1;
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            B_(i,j)=phi0(coefficient((sBXProm_(0,i)**sBYProm_(0,j))_(0,0), bezDetRed));
        );
    );
    return matrix(B);
);






