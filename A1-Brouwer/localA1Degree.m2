-- Given an endomorphism of affine space f=(f1,...,fn) given as a list of polynomials, return the Bezoutian corresponding to the endomorphism
--for notation: input domain is kk[x_1..x_n] while Bezoutian computed in kk[X_1..Y_n]

path = append(path, "/home/macaulay/A1-Brouwer/");
load "GW-type.m2"
load "localAlgebraBasis.m2"

debugging = false;


localA1Degree = method()
localA1Degree (List, Ideal) := (Matrix) => (Endo,p) -> (
    -- Endo is the list {f_1, f_2, ..., f_n} of polynomials kk^n -> kk^n
    
    -- n is the number of polynomials
    n := #Endo; -- n is the no. of polynomials
    

    -- Get the underlying field, assert it is a field
    kk := coefficientRing(ring(Endo#0)); 
    
       
    if isField(kk) == false then(
	kk = toField(kk);
	);
    
    if debugging==true then(
    	print("kk is " | toString(kk));
    	print("kk is a field? " | toString(isField(kk)));
	);
    
    -- Let S = k[x_1..x_n] be the ambient polynomial ring
    S:=ring(Endo#0);
    
    
    -- Create the ideal J.  It has prop that  k[x_1 .. x_n]_p/I_p is isom to k[x_1 .. x_n]/J    
    J:=(ideal Endo):saturate(ideal Endo,p);
    
    -- Get the dimension of the local algebra Q_p(f)
    localFormRank := numColumns(basis(S/J));
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then error "Error: morphism does not have isolated zeroes";
    
    -- Check the number of variables matches the number of polynomials
    if not #(gens S) == n then error "Error: the number of variables does not match the number of polynomials.";
    
    -- If the field is CC, just output gwClass of an identity matrix of rank = localFormRank
    if (kk === CC or instance(kk,ComplexField)) then(
	return gwClass(matrix(mutableIdentity(CC,localFormRank)))
    );
    
    -- If the field is RR, ask the user to run it over QQ instead, then simplify over RR
    if (kk === RR or instance(kk,RealField)) then error "Error: globalA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output a GrothendieckWittClass. Then extract the matrix, base change it to RR, and then run simplifyForm().";
    
    
    -- Create internal rings/matrices
    
    -- Initialize a polynomial ring in X_i's and Y_i's to compute the Bezoutian in
    R := kk[X_1..Y_n];
    
    
    -- Create an (n x n) matrix D which will be populated by \Delta_{ij} in the paper    
    D := "";
    try D = mutableMatrix id_((frac R)^n) else D= mutableMatrix id_(R^n);
    
    
    if debugging==true then(
    	print("R is " | toString(R));
    	print("D is");
    	print(D);
	);    
    
    
    -- Iterate through the entries of the matrix D and populate it with the Delta_{ij}'s
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(-- TODO let's rewrite this to just run 1 to n instead (?)
	    
	    if debugging == true then(
	    	print("i=" | toString(i) | ", j="  | toString(j));
		);
	    
	    
	     -- Creates the list {x_1 => Y_1, ..., x_j-1 => Y_j-1, x_j => X_j, ..., x_n => X_n}
	    targetList1 := toList (Y_1..Y_j|X_(j+1)..X_n);
	    
	     -- Creates the list {x_1 => Y_1, ..., x_j => Y_j, x_j+1 => X_j+1, ..., x_n => X_n}
	    targetList2 := toList (Y_1..Y_(j+1)| X_(j+2)..X_n);
            
	    -- Build the quantity f_i(Y_1,..X_{j+1},..X_n) - f_i(Y_1..Y_{j+1},..,X_n)
	    numeratorD = ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i));
	    
	    -- Divide by X_{j+1}-Y_{j+1}	    
	    D_(i,j)= numeratorD/(X_(j+1)-Y_(j+1)); -- this is \Delta_{i,j} from the paper
	); 
    );

    -- TESTING
     if debugging==true then(
    	 print("\nfracFieldBezDet is " | toString(det(D)));
    	 print("\nRing of fracFieldBezDet is " | toString(ring det(D)));
    	 print("\nnumerator is " | toString(numerator(det(D))));
    	 print("\ndenominator is " | toString(denominator(det(D))));
	 print("\nR is " | toString(R));
    	 print("\nlift is " | toString(numerator(det(D)), R));
    	 print("lift denom to coeff ring is " | toString(lift(denominator(det(D)),coefficientRing R)));
	 );
    
    -- Set up the local variables bezDet and bezDetR
    bezDet:="";
    bezDetR:="";
    
    
    -- The determinant of D is interpreted as living in Frac(k[x_1..x_n]),
    -- so we can try to lift it to k[x_1..x_n]           
    if liftable(det(D),R) == true then(
	bezDetR = lift(det(D),R);
	);
    
    -- In some computations, applying lift(-,R) won't work. So we instead lift
    -- the numerator and then divide by a lift of the denominator (which will be
    -- a scalar) to the coefficient ring k
    if not liftable(det(D),R) == true then(
	bezDet = lift(numerator(det(D)), R) / lift(denominator(det(D)),coefficientRing R);
    	bezDetR = lift(bezDet, R);
	);    
    
    -- Define formal variables X_i, Y_i that replace x_i
    RX:=kk[X_1..X_n]; 
    RY:=kk[Y_1..Y_n];
    
    -- mapxtoX replaces all instances of x_i with X_i. mapxtoY does the same but with Y_i's
    mapxtoX:= (map(RX,S,toList(X_1..X_n)));
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n)));

    -- Find standard basis and define local quotient ring
    list1 := (apply(toList(0..n-1), i-> mapxtoX(Endo_i))); -- list (f_1(X), ..., f_n(X))
    list2 := (apply(toList(0..n-1), i-> mapxtoY(Endo_i))); -- list (f_1(Y), ..., f_n(Y))
    
    -- Apply localAlgebraBasis to compute standard basis for localization
    standBasisX :=localAlgebraBasis(list1,mapxtoX p);
    standBasisY := localAlgebraBasis(list2,mapxtoY p); 
    

    
    -- Take the ideal sum of J in the X_i's with J in the Y_i's
    localIdeal :=sub(mapxtoX(J),R)+sub( mapxtoY(J),R);
    
    
    
    
    -- Here we're using that (R/I) \otimes_R (R/J) = R/(I+J) in order to express Q_p(f) \otimes Q_p(f), where X's are the variables in first term, Y's are variables in second par
    Rquot:=R/localIdeal;

    if debugging==true then(
    	print("\nring of localIdeal is " | toString(ring localIdeal));
    	print("\nRquot is " | toString(Rquot));
	);
    
    
    -- Move the standard bases to the quotient ring
    sBXProm :=apply(toList(0..#standBasisX-1),i-> sub(standBasisX_i,Rquot));
    sBYProm :=apply(toList(0..#standBasisY-1),i-> sub(standBasisY_i,Rquot));
   
    -- Now reduce the bezDet determinant subject ot the local ideal in both the X's, Y's.
    bezDetRed := bezDetR % localIdeal;
    
    if debugging==true then(
    	print("\nbezdetred is " | toString(bezDetRed));
	);
    
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
