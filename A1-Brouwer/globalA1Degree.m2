--Nikita, Gabriel, Jordy
-- Given an endomorphism of affine space f=(f1,...,fn) given as a list of polynomials, return the Bezoutian corresponding to the endomorphism
--for notation: input domain is kk[x_1..x_n] while Bezoutian computed in kk[X_1..Y_n]
load "./GW-type.m2"
load "./rankGlobalAlgebra.m2"

globalA1DegreeCC = method()
globalA1DegreeCC (List) := (GrothendieckWittClass) => (Endo) ->(
    print("globalA1DegreeCC CALLED");
    if dim ideal(Endo) > 0  then error "Error: morphism does not have isolated zeroes";
    S:=ring(Endo#0);
    rankAlgebra:=rankGlobalAlgebra(Endo);
    return gwClass(matrix(mutableIdentity(CC,rankAlgebra)));  
    );

globalA1Degree = method()
globalA1Degree (List) := (GrothendieckWittClass) => (Endo) -> (
    -- Endo is the list {f_1, f_2, ..., f_n} of polynomials kk^n -> kk^n
    
    -- n is the number of polynomials
    n := #Endo;
    
    -- Get the underlying field    
    kk := coefficientRing(ring(Endo#0));
    if isField(kk) == false then(
    	kk = toField(kk);
    	);
    
    -- Let S = k[x_1..x_n] be the ambient polynomial ring
    S:=ring(Endo#0);    
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then error "Error: morphism does not have isolated zeroes";
    
    -- Check the number of variables matches the number of polynomials
    if not #(gens S) == n then error "Error: the number of variables does not match the number of polynomials.";
    
    -- If the field is CC, run the much quicker globalA1DegreeCC method
    if (kk === CC or instance(kk,ComplexField)) then(
	print("Got here");
    return globalA1DegreeCC(Endo)
    );
    
    -- If the field is RR, ask the user to run it over QQ instead, then simplify over RR
    if (kk === RR or instance(kk,RealField)) then error "Error: globalA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output a GrothendieckWittClass. Then extract the matrix, base change it to RR, and then run simplifyForm().";
    
    
    -- Create internal rings/matrices
    
    -- Initialize a polynomial ring in X_i's and Y_i's to compute the Bezoutian in
    R := kk[X_1..Y_n];
    
    -- Create an (n x n) matrix D which will be populated by \Delta_{ij} in the pape
    D := "";
    try D = mutableMatrix id_((frac R)^n) else D= mutableMatrix id_(R^n);
    
    
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(
	    -- iterate through the entries of the matrix D and populate it with the following information ...
        -- create the list {x_1 => Y_1, ..., x_j-1 => Y_j-1, x_j => X_j, ..., x_n => X_n}
	    targetList1 := toList (Y_1..Y_j|X_(j+1)..X_n);
        -- create the list {x_1 => Y_1, ..., x_j => Y_j, x_j+1 => X_j+1, ..., x_n => X_n}
	    targetList2 := toList (Y_1..Y_(j+1)| X_(j+2)..X_n); 
        -- map f_i(x) to f_i(Y_1, ..., Y_j-1, X_j, ..., X_n) resp.
        numeratorD = ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i)); 
        -- this is \Delta_{i,j} from the paper 
	    D_(i,j)= numeratorD/(X_(j+1)-Y_(j+1)); 
	    ); 
        );


    -- The determinant of D is interpreted as living in Frac(k[x_1..x_n]).
    -- Applying lift(-,R) won't work here, so we lift the numerator and then
    -- divide out by a lift of the denominator (which will be a scalar) to the
    -- coefficient ring 
    fracFieldBezDet := det(D);
    bezDet := lift(numerator(det(D)), R) / lift(denominator(det(D)),coefficientRing R);
    bezDetR := lift(bezDet, R);
    
    RX:=kk[X_1..X_n]; 
    RY:=kk[Y_1..Y_n];
    -- defines the map f_i(x_1, ..., x_n) to f_i(X_1, ..., X_n)
    mapxtoX:= (map(RX,S,toList(X_1..X_n))); 
    -- defines the map f_i(x_1, ..., x_n) to f_i(Y_1, ..., Y_n)
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n))); 
    -- Compute the standard basis of kk[X_1, ..., X_n]/(f_1, ..., f_n)
    standBasisX := basis (RX/(ideal (leadTerm (mapxtoX ideal Endo)))); 
    
    -- Compute the standard basis of kk[Y_1, ..., Y_n]/(f_1, ..., f_n)
    standBasisY := basis (RY/(ideal (leadTerm (mapxtoY ideal Endo)))); 
    -- defines an ideal (f_1(X), ..., f_n(X))
    id1 := (ideal apply(toList(0..n-1), i-> mapxtoX(Endo_i))); 
    -- defines an ideal (f_1(Y), ..., f_n(Y))
    id2 := (ideal apply(toList(0..n-1), i-> mapxtoY(Endo_i))); 
    -- takes the sum of the ideals in the ring kk[X_1..Y_n]
    promotedEndo :=sub(id1,R)+sub(id2,R); 
    -- defines the quotient ring kk[X_1..Y_n]/(f(X),f(Y)) which is Q(f)\otimes_{k}Q(f) in the paper
    Rquot:= R/promotedEndo; 
    -- moves the standard bases to the quotient ring
    sBXProm :=sub(standBasisX,Rquot); 
    -- moves the standard bases to the quotient ring 
    sBYProm :=sub(standBasisY,Rquot); 
    
    --reduces bezDet mod promotedEndo
    bezDetRed := bezDetR % promotedEndo; --reduces bezDet mod promotedEndo 

    
    ------------------------
    -- define a ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring
    phi0 := map(kk,Rquot,(toList ((2*n):0))); 
    m:= numColumns(sBXProm);
    B:= mutableMatrix id_(kk^m);
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            B_(i,j)=phi0(coefficient((sBXProm_(0,i)**sBYProm_(0,j))_(0,0), bezDetRed));
        );
    );
    return matrix(B);

);
