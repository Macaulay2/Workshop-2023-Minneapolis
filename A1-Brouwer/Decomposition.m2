

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



