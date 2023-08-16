-- Given a form q over QQ of anisotropic dimension d>=4, returns a form <a> so that q+<-a> has anisotropic dimension d-1
-- This is Koprowski/Rothkegel's Algorithm 5 in the case of QQ
QQanisotropicDimension4 = method()
QQanisotropicDimension4 (GrothendieckWittClass) := (GrothendieckWittClass) => beta ->(
    if not (anisotropicDimensionQQ(beta) >= 4) then error "anisotropic dimension of inputted form is not >=4";
    
    -- If the signature is non-negative then return <1>
    if signature(beta) >= 0 then(
	return gwClass(matrix(QQ,{{1}}))
	);
    
    -- Otherwise return <-1>
    if signature(beta) < 0 then(
	return gwClass(matrix(QQ,{{-1}}))	        
        );
	
    );

-- Given a form q over QQ of anisotropic dimension 3, returns a form <a> so that q+<-a> has anisotropic dimension 2
-- This is Koprowski/Rothkegel's Algorithm 7 in the case of QQ
QQanisotropicDimension3 = method()
QQanisotropicDimension3 (GrothendieckWittClass) := (GrothendieckWittClass) => beta ->(
    d := integralDiscriminant(beta);
    
    -- Build lists of primes where the p-adic valuation of the discriminant is even or is odd
    L1:={};
    L2:={};
    S1:={};
    S2:={};
    for p in relevantPrimes(beta) do(
	if odd PadicValuation(d,p) then(
	    L1 = append(L1,p);
	    S1 = append(S1,d-1);
	    );
	if even PadicValuation(d,p) then(
	    L2 = append(L2,p^2);
	    S2 = append(L2,p)
	    );
	);
    
    -- We are looking for an element which is equivalent to d-1 mod p for each p in L1, and equivalent to p mod p^2 for each p in L2
    -- we use the chineseRemainder method from the "Parametrization" package to find such an element
    alpha := chineseRemainder(S1 | S2, L1 |L2);
    a := squarefreePart(alpha);
    
    return gwClass(QQ,{{a}});

    );

-- Constructs the anisotropic part of a form with anisotropic dimension 2
-- 
-- QQanisotropicDimension2 = method()
-- QQanisotropicDimension2 (GrothendieckWittClass) := (GrothendieckWittClass) => beta ->(
--     -- If the witt Index isn't 0 mod 4, add on some hyperbolic forms so that 
--     n := numRows beta.matrix;
--     wittIndexBeta := n - anisotropicDimensionQQ(beta);
--     q:= beta;
--     -- If the Witt Index isn't 0 mod 4 it must be 2 mod 4, so we add on a hyperbolic form
--     if (not wittIndexBeta % 4 == 0) then(
-- 	q = gwAdd(beta, hyperbolicForm(QQ));
-- 	);
    
--     d:= integralDiscriminant(q);
    
--     W:= matrix(QQ,{{}});
--     if (d <= 0) then(
-- 	W = W | matrix(QQ,{{1/2 - signature(q)/4}});
-- 	);
    
--     for p in relevantPrimes(beta) do(
-- 	W = W | matrix(QQ,{{(1 - (HasseWittInvariant(q,p)))/2}});	
-- 	);    
--    );






WittDecomp = method()
WittDecomp (Matrix) := (ZZ,Matrix) => (A) -> (
    k := ring A;   
  
    -- Add error in case the base field is RR or CC
    if (instance(k,InexactFieldFamily) or instance(k,RealField) or instance(k,ComplexField)) then error "Error: base field is inexact, use WittDecompInexact() instead";
    
    -- Rank of matrix
    n := numRows(A);
    x := symbol x;
    -- Local variable ring
    R := k[x_0..x_(n - 1)];
    
    -- Create quadratic from f from matrix A
    f := sum (
        for i from 0 to (n - 1) list (
            sum (for j from 0 to (n - 1) list (A_(i,j)*x_i*x_j))
            )
        );
 
    -- will use rationalPoints package to see if f has a zero.  If so, A is isotropic.  
    -- rationalPoints package seeks a zero upto variable bound
    
    use k;
    solnPt := new MutableList;
    solnFound := false;
    for bound from 1 to 10 do (
        solns := new MutableList from rationalPoints(ideal(f),Bound => bound);
        if (#solns > 1) then(
            if solns#0 == toList(n:0) then (solnPt = solns#1) else (solnPt = solns#0);
            solnFound = true;
            break;
        );
    );


    -- If no solutions found, then we assume the form is anisotropic
    if (not solnFound) then (return (0,A));
    
    -- If solution is found for a rank 2 form, then the form is purely hyperbolic
    if ((n == 2) and  (det(A) == 0))then (error "Matrix singular, run WittDecompGeneral instead" );
    if (n == 2) then (return (1,matrix(k,{{}})));

    -- If found a solution, record it as row matrix z. Then find vector y such that < z,y > <>0.  
    -- z and y then generate a hyberbolic plane. 
    -- z as a row matrix
    z := matrix{toList solnPt};
    zA := z*A;
    y := new MutableMatrix from matrix{toList(n:(0/1))};
    for i from 0 to (n - 1) do (
        if (zA_(0,i) != 0) then (y_(0,i) = 1; break;);
    );
    
    -- Now z and y span a copy of |H in the bilinear form
    -- We need to find a basis of vectors orthogonal (wrt bilinear form) to x and y 
    -- An (n - 2) x n matrix.
    orthoComp := gens kernel((z||matrix(y))*A);
    
    -- Now recursively apply WittDecomp to orthoComp^T*A*orthoComp a (n - 2) x (n - 2) Gram matrix
    subComputation := WittDecomp(transpose(orthoComp)*A*orthoComp);
    
    -- subComputation_0 gives number of hyperbolic forms in (n - 2) x (n - 2) subform  
    -- 1+ subComputation_0 is the number of hyperbolic forms in A
    -- subComputation_1 is the anisotropic part of A and (also) the subform part
    
    return (1+subComputation_0, subComputation_1);
    )

---------------------------------------
-- Simplifying a form
---------------------------------------

sumDecompositionVerbose = method()
sumDecompositionVerbose (GrothendieckWittClass) := (GrothendieckWittClass, String) => beta -> (
    -- Check if the diagonalForm has already been computed, if so recall it from the cache    
    gamma := diagonalForm(beta);
    
    -- Get the matrix, base field, and rank of the form
    A := gamma.matrix;
    k := baseField(gamma);
    n := numRows(gamma.matrix);
    
    -- If the form has no rank, stop here and return the zero form. This gets rid of fringe cases later on
    if n == 0 then return (gamma,"zero form");
    
    -- Define the hyperbolic form over our field
    H := matrix(k,{{1,0},{0,-1}});
        
    -----------------------
    -- Field by field cases
    -----------------------
        
    --In the case of C:
    if (k === CC or instance(k,ComplexField)) then(
        if n == 2 then(
            return (H, "H");
        );
        if even n and n > 2 then(
        
            -- The form is (n/2)H
            numHypFormsCCEven := n/2;
        
            simplifiedFormCCEven := H;
        
            for i in 1..sub((numHypFormsCCEven-1),ZZ) do(
                simplifiedFormCCEven = simplifiedFormCCEven ++H;
		);
        
            outputStringCCEven := toString(numHypFormsCCEven) | "H";
            return (gwClass(simplifiedFormCCEven),outputStringCCEven)
            );
        if odd n then(
        
            -- The form is floor(n/2)H + <1>
            numHypFormsCCOdd := floor(n/2);
            simplifiedFormCCOdd := H;
        
            for i in 1..sub((numHypFormsCCOdd-1),ZZ) do(
                simplifiedFormCCOdd = simplifiedFormCCOdd ++H;
		);
        
            simplifiedFormCCOdd = simplifiedFormCCOdd ++ matrix(k,{{1}});
        
            outputStringCCOdd := toString(numHypFormsCCOdd) | "H + <1>";
            return (gwClass(simplifiedFormCCOdd),outputStringCCOdd)
            );
	);
	-- End CC case
    
    -- If the field is R, look at the sign of elements along the diagonal
    if (k === RR or instance(k,RealField)) then(
        posEntries := 0; --for loop counts the number of positive diagonal entries of diagA
        negEntries := 0; --for loop counts the number of negative diagonal entries
	for i from 0 to (n-1) do(
            if A_(i,i)>0 then(
                posEntries=posEntries+1;
            	);
            if A_(i,i)<0 then(
                negEntries=negEntries+1;
            	);
            );
            
	-- Number of hyperbolic forms is the number of positive and negative entries
        wittIndexRR := min(posEntries,negEntries);
    	
	-- Make an empty matrix and string to add output to
	simplifiedFormRR := matrix(k,{{}});
	outputStringRR := "";
	
	-- Add hyperbolic forms to output
	for i in 1..(wittIndexRR) do(
	    simplifiedFormRR = safeBlockSum(simplifiedFormRR,H);
            );
	
    if wittIndexRR == 1 then(
        outputStringRR = outputStringRR | "H";
    );

	if wittIndexRR > 1 then(
	    outputStringRR = outputStringRR | toString(wittIndexRR) | "H";
	    );
	
	-- Look at what's left over
	posEntries = posEntries - wittIndexRR;
	negEntries = negEntries - wittIndexRR;
	
	-- Build the anisotropic part
	anisotropicPart := safeBlockSum(matrix(mutableIdentity(k,posEntries)),((-1)*matrix(mutableIdentity(k,negEntries))));
	
	-- Add on (safely) the anisotropic part
	simplifiedFormRR = safeBlockSum(simplifiedFormRR,anisotropicPart);
	
	-- Amend the string output with anisotropic info
	if posEntries == 1 then(
	    outputStringRR = outputStringRR | " + " | "<1>";
	    );
	if posEntries > 1 then(
        outputStringRR = outputStringRR | " + " | toString(posEntries) | "<1>";
    );
    if negEntries == 1 then(
        outputStringRR = outputStringRR | " + " | "<-1>";
    );
	if negEntries > 1 then(
	    outputStringRR = outputStringRR | " + " | toString(negEntries) | "<-1>";
	    );
	
	return (gwClass(simplifiedFormRR), outputStringRR);
	
	);
	-- End RR case
	-- TODO are we getting the string "4H" or whatever in the output string?
       
       
    -- Finite field case
    
    if instance(k, GaloisField) then(
	
	-- Take the diagonal entries to a list
	diagEntries := {};
	for i from 0 to sub((n-1),ZZ) do(
	    diagEntries = append(diagEntries, A_(i,i));
	    );
	
	-- Start counting squares and nonsquares
        numSquares:= 0;
	numNonSquares:= 0;
	
	-- Define a new variable which we will set equal to a nonsquare representative
	nonSquareRepresentative := sub(0,k);
	
	
	-- Check if each elt is a square or not
	for x in diagEntries do(
	    
	    -- If it is a square
	    if legendreBoolean(x) then(
		numSquares = numSquares + 1;
		);
	    
	    -- If it is not a square
	    if not legendreBoolean(x) then(
		--Pick x to be a representative of the non-square class in F_q^x / (F_q^x)^2
		nonSquareRepresentative = x;
		numNonSquares = numNonSquares + 1;
		);
	    );
	
	
	-- If -1 is a square class in k, then 2<1> = H
	if legendreBoolean(sub(-1,k)) then(
	    
	    -- Set up output matrix and string
	    simplifiedFormGFSquare := matrix(k,{{}});
	    outputStringGFSquare := "";	    
	    
	    -- Number of hyperbolic forms
            wittIndexGFSquare := floor(numSquares/2) + floor(numNonSquares/2);
	    
	    -- If there are hyperbolic forms add them to the matrix and string
        if wittIndexGFSquare == 1 then(
        outputStringGFSquare = outputStringGFSquare | "H";
        simplifiedFormGFSquare = H;
        );

	    if wittIndexGFSquare > 1 then(
		outputStringGFSquare = outputStringGFSquare | toString(wittIndexGFSquare) | "H";
		for i in 1..(wittIndexGFSquare) do(
		    simplifiedFormGFSquare = safeBlockSum(simplifiedFormGFSquare,H);
		    );
		);
	    
	    -- There will be at most one extra square or one extra nonsquare left over	
	    if odd numSquares then(
		outputStringGFSquare = outputStringGFSquare | " + <1>";		
		simplifiedFormGFSquare = safeBlockSum(simplifiedFormGFSquare ,matrix(k,{{1}}));		
		);
	    
	    if odd numNonSquares then(
		outputStringGFSquare = outputStringGFSquare | " + <" | toString(nonSquareRepresentative) | ">";		
		simplifiedFormGFSquare = safeBlockSum(simplifiedFormGFSquare,((-1)*matrix(k,{{nonSquareRepresentative}})));
		
		);
		
	    return(gwClass(simplifiedFormGFSquare),outputStringGFSquare)
	    
	    
	    ); -- End "-1 is a square in GF"
	
	-- If -1 is not a square then it is a representative for the nonsquare class	  	   
	   if not legendreBoolean(sub(-1,k)) then (
	      	            
	       -- Set up output matrix and string
	       simplifiedFormGFNonSquare := matrix(k,{{}});
	       outputStringGFNonSquare := "";
	       
	       
	       -- Number of hyperbolic forms
	       wittIndexGFNonSquare := min(numSquares,numNonSquares);	       
	       
	       -- Look at any remaining squares or nonsquares. One of these must be zero
	       numSquares = numSquares - wittIndexGFNonSquare;
	       numNonSquares = numNonSquares - wittIndexGFNonSquare;
	       
	       
	       -- Add on hyperbolic part
           if wittIndexGFNonSquare == 1 then(
           simplifiedFormGFNonSquare = H;
           outputStringGFNonSquare = outputStringGFNonSquare | "H";
           );

	       if wittIndexGFNonSquare > 1 then(
		   for i in 1..(wittIndexGFNonSquare) do(
		       simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, H);
		       
		       );
		   outputStringGFNonSquare = outputStringGFNonSquare | toString(wittIndexGFNonSquare) | "H";
		   
		   );
	       
	       -- Add any anisotropic part
           if numSquares == 1 then(
           simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, matrix(mutableIdentity(k,numSquares)));
		   outputStringGFNonSquare = outputStringGFNonSquare | " + " | "<1>";
           );
	       if numSquares > 1 then(
		   simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, matrix(mutableIdentity(k,numSquares)));
		   outputStringGFNonSquare = outputStringGFNonSquare | " + " | toString(numSquares) | "<1>";
		   );
		   if numNonSquares == 1 then(
		   simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, ((-1)*matrix(mutableIdentity(k,numNonSquares))));
		   outputStringGFNonSquare = outputStringGFNonSquare | " + " |  "<-1>";
		   );
	       if numNonSquares > 1 then(
		   simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, ((-1)*matrix(mutableIdentity(k,numNonSquares))));
		   outputStringGFNonSquare = outputStringGFNonSquare | " + " | toString(numNonSquares) | "<-1>";
		   );
	       
	       return (gwClass(simplifiedFormGFNonSquare), outputStringGFNonSquare)		   
	       
	       ); -- End "-1 is not a square in GF"
	
	    
	  
			
      
	   
	   
	   );-- End finite field case
       
          -- Number field case
    
    -- NOTE: I wanted to do member(QQ,k.baseRings) but this returns false.
    -- 	   This is because baseField(..) will return a ring
    
    
    -- We're making the number field case a catch-all case now
--    if member(QQ, k.baseRings) and (isField(k)) then(
	-- Set up output form and matrix
	simplifiedFormQQ := matrix(k,{{}});
        outputStringQQ := "";
	
	-- Get number of confirmed hyperbolic forms and remainder from WittDecomp
	(numHypForms,B) := WittDecomp(A);
	
	-- Add any hyperbolic forms if they exist
    if numHypForms = 1 then(
        simplifiedFormQQ = H;
        outputStringQQ = outputStringQQ | "H + potentially anisotropic part";
    );
	if numHypForms > 1 then(
	    for i in 1..(numHypForms) do(
		
		simplifiedFormQQ = safeBlockSum(simplifiedFormQQ, H);
		outputStringQQ = outputStringQQ | toString(numHypForms) | "H + potentially anisotropic part";
		
		);
	    );
	
	if numHypForms == 0 then(
	    outputStringQQ = outputStringQQ | " potentially anisotropic part";
	    
	    );
	 simplifiedFormQQ = safeBlockSum(simplifiedFormQQ, B);
	 
	 return (gwClass(simplifiedFormQQ),outputStringQQ)
		
--	);-- End number field case   
    );


sumDecomposition = method()
sumDecomposition (GrothendieckWittClass) := (GrothendieckWittClass) => beta -> (
    beta.cache.diagonalForm = (sumDecompositionVerbose(beta))_0;
    return (sumDecompositionVerbose(beta))_0
);

sumDecompositionString = method()
sumDecompositionString (GrothendieckWittClass) := (String) => beta -> (
    return (sumDecompositionVerbose(beta))_1
);



