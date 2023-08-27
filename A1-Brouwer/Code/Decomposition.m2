
-- Input: A form 1 over QQ of anisotropic dimension d >= 4
-- Output: A form < a > so that q + < a > has anisotropic dimension d - 1

-- Note: This is Koprowski/Rothkegel's Algorithm 5 in the case of QQ

qQanisotropicDimension4 = method()
qQanisotropicDimension4 (GrothendieckWittClass) := (GrothendieckWittClass) => beta ->(
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

-- Input: A form q over QQ of anisotropic dimension 3
-- Output: A form < a > so that q + < -a > has anisotropic dimension 2

-- Note: This is Koprowski/Rothkegel's Algorithm 7 in the case of QQ

qQanisotropicDimension3 = method()
qQanisotropicDimension3 (GrothendieckWittClass) := (GrothendieckWittClass) => beta ->(
    d := integralDiscriminant(beta);
    
    -- Build lists of primes where the p-adic valuation of the discriminant is even or is odd
    L1 := {};
    L2 := {};
    S1 := {};
    S2 := {};
    for p in relevantPrimes(beta) do(
	if odd padicValuation(d,p) then(
	    L1 = append(L1,p);
	    S1 = append(S1,d-1);
	    );
	if even padicValuation(d,p) then(
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
qQanisotropicDimension2 = method()
qQanisotropicDimension2 (GrothendieckWittClass) := (GrothendieckWittClass) => beta ->(
    n := numRows beta.matrix;

    -- Shortcut: if the form has anisotropic dimension 2 and the form is dimension 2, return the form itself
    -- if (n==2) then(
    -- 	return beta;
    -- 	);
        
    
    -- Step 1: We want the Witt index to be 0 mod 4 in their terminology --- note they define the Witt index to be
    -- the integer w so that q = wH + q_a. This is not the same as the dimension of a maximal totally isotropic subspace
    w := (n - anisotropicDimensionQQ(beta))/2;
    w = sub(w,ZZ);
    q := beta;
    if ((w % 4) != 0) then(
	w = w % 4;
	q = gwAdd(q, hyperbolicForm(QQ,2*(4-w)));
	n = n + 2*(4-w);
	);
    -- Step 2: Compute discriminant (note they use a signed version of the discriminant in their algorithm)
    d:= ((-1)^(n*(n-1)/2))*integralDiscriminant(q);
    
    -- Step 3: Take relevant primes plus dyadic ones
    L := relevantPrimes(beta);
    if not member(2,L) then(
	L = append(L,2);
	);
    
    -- Start the loop at p=2
    p:=2;
    solnFound := false;
    
    
    while not solnFound do(
	r := #L;
    	-- Step 5c: Make a vector of exponents of Hasse invariants
	W := mutableMatrix(QQ,r,1);
	for i from 0 to (r-1) do(
	    W_(i,0) = (1 - (hasseWittInvariant(q,L_i)))/2;
	    );
       	

	-- Step 5b: 
	W = matrix(W);
    	if (d < 0) then(
	    if not (abs(signature(q)) == 2) then (error "signature isn't pm 2");
	    
	    if (signature(q) == 2) then (
		W = matrix(QQ,{{0}}) || W;
		);
	    if (signature(q) == -2) then(
		W = matrix(QQ,{{1}}) || W;
	    );
	);
        
    	-- Step 5e: Make a matrix of Hilbert symbols
    	B := mutableIdentity(QQ,r);
    	for i from 0 to (r-1) do(
	    for j from 0 to (r-1) do(
	    	B_(i,j) = (1 - hilbertSymbol(L_j, d, L_i))/2;
	    	);
	    );
	B = matrix(B);
    	
	-- Step 5d: Append a zero column on the front if the discriminant is negative
    	if (d < 0) then(
	    zeroVec := mutableMatrix(QQ,1,r);
	    for i from 0 to (r-1) do(
	    	zeroVec_(0,i) = 0
	    	);
	    B = matrix(zeroVec) || B;
	    );
        kk := GF(2);
    	W = matrix(kk,entries(W));
    	B = matrix(kk,entries(B));

	
	if (class(solve(B,W)) === Matrix) then(
	    X := solve(B,W);
	    solnFound = true;
	    break;
	    )
	else(
	    p = nextPrime(p+1);
	    while (member(p,L)==true) do(
		p = nextPrime(p+1);
		);

	    L = append(L,p);
	    );
	);
  
    alpha := sub(1,ZZ);
    for j from 0 to (r-1) do(
	alpha = alpha * ((L_j)^(sub(X_(j,0),ZZ)));
	);
    return diagonalClass(QQ,(alpha, -alpha*d))
    
   );



-- Input: Any form over QQ
-- Output: Its anisotropic part
qQanisotropicPart = method()
qQanisotropicPart (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    beta = integralDiagonalRep(beta);
    
    n := numRows(beta.matrix);
    d := anisotropicDimension(beta);
    
    -- If the form is anisotropic 
    if n == d then(return beta);
    
    -- Initialize an empty quadratic form
    outputForm := diagonalClass(QQ,());
    alpha := 1;
    
    
    while d>=4 do(
	d = anisotropicDimension(beta);
	outputForm = gwAdd(outputForm,qQanisotropicDimension4(beta));
	alpha = (qQanisotropicDimension4(beta).matrix)_(0,0);
	
	beta = gwAdd(beta, diagonalClass(QQ,((-1)*alpha)));
	
	);
    
    if d==3 then(
	outputForm = gwAdd(outputForm,qQanisotropicDimension3(beta));
	alpha = (qQanisotropicDimension3(beta).matrix)_(0,0);
	
	beta = gwAdd(beta, diagonalClass(QQ,((-1)*alpha)));
	
	);
    
    if d==2 then(
       outputForm = gwAdd(outputForm, qQanisotropicDimension2(beta));
       );
    
    if d==1 then(
	outputForm = gwAdd(outputForm, diagonalClass(QQ,(integralDiscriminant(beta))));
	);
    
    return outputForm;
    
    
    );

-- Input: A symmetric matrix representing a quadratic form or a GrothendieckWittClass; over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: A symmetric matrix or GrothendieckWittClass that is the anisotropic part of the input

anisotropicPart = method()
anisotropicPart (Matrix) := (Matrix) => (A) -> (
    k := ring A;
    -- Ensure base field is supported
    if not (k === CC or instance(k,ComplexField) or k === RR or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and k.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    -- Ensure underlying matrix is symmetric
    if (transpose(A) != A) then (
        error "Underlying matrix is not symmetric";
	);
    diagA := congruenceDiagonalize(A);
    -- Over CC, the anisotropic part is either the rank 0 form or the rank 1 form, depending on the anisotropic dimension
    if (k === CC or instance(k,ComplexField)) then (
        if (anisotropicDimension(A)==0) then (
            return (diagonalMatrix(CC,{}));
            )
        else (
            return (matrix(CC,{{1}}));
            );
        )
    --Over RR, the anisotropic part consists of the positive entries in excess of the number of negative entries, or vice versa
    else if (k === RR or instance(k,RealField)) then (
        posEntries := numPosDiagEntries(diagA);
        negEntries := numNegDiagEntries(diagA);
        if (posEntries > negEntries) then (
            return (id_(RR^(posEntries-negEntries)));
            )
        else if (posEntries < negEntries) then (
            return (-id_(RR^(negEntries-posEntries)));
            )
        else (
            return (diagonalMatrix(RR,{}));
            );
        )
    -- Over QQ, call anisotropicPartQQ
    else if (k === QQ) then (
        return (qQanisotropicPart(gwClass(nondegeneratePart(A)))).matrix;
        )
    -- Over a finite field, if the anisotropic dimension is 1, then the form is either <1> or <e>, where e is any nonsquare representative, and if the anisotropic dimension is 2 then the form is <1,-e>
    else if (instance(k, GaloisField) and k.char != 2) then (
        if (anisotropicDimension(A)==1) then (
            return (matrix(k,{{sub((-1)^((numNonzeroDiagEntries(diagA)-1)/2),k)*det(nondegeneratePartDiagonal(diagA))}}));
            )
        else if (anisotropicDimension(A)==0) then (
            return (diagonalMatrix(k,{}));
            )
        else (
            return (matrix(k,{{1,0},{0,sub((-1)^((numNonzeroDiagEntries(diagA)-2)/2),k)*det(nondegeneratePartDiagonal(diagA))}}));
            );
        )
    -- We should never get here
    else error "Problem with base field"
    )


anisotropicPart (GrothendieckWittClass) := (GrothendieckWittClass) => (alpha) -> (
    return (gwClass(anisotropicPart(alpha.matrix)));
    )


-- Input: A matrix whose base field is not inexact
-- Output: The number of hyperbolic forms in A, and the anisotropic part of A and (also) the subform part

wittDecomp = method()
wittDecomp (Matrix) := (ZZ,Matrix) => (A) -> (
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
    if n>1 then(
    for bound from 1 to 10 do (
        solns := new MutableList from rationalPoints(ideal(f),Bound => bound);
        if (#solns > 1) then(
            if solns#0 == toList(n:0) then (solnPt = solns#1) else (solnPt = solns#0);
            solnFound = true;
            break;
        );
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
    subComputation := wittDecomp(transpose(orthoComp)*A*orthoComp);
    
    -- subComputation_0 gives number of hyperbolic forms in (n - 2) x (n - 2) subform  
    -- 1+ subComputation_0 is the number of hyperbolic forms in A
    -- subComputation_1 is the anisotropic part of A and (also) the subform part
    
    return (1+subComputation_0, subComputation_1);
    )

---------------------------------------
-- Simplifying a form
---------------------------------------

-- Input: A Grothendieck-Witt class beta over a field k
-- Output: A simplified diagonal representative of beta

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
	(numHypForms,B) := wittDecomp(A);
	
	-- Add any hyperbolic forms if they exist
    if numHypForms == 1 then(
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

-- Input: A Grothendieck-Witt class beta over a field k
-- Output: A simplified diagonal representative of beta

sumDecomposition = method()
sumDecomposition (GrothendieckWittClass) := (GrothendieckWittClass) => beta -> (
    beta.cache.diagonalForm = (sumDecompositionVerbose(beta))_0;
    return (sumDecompositionVerbose(beta))_0
);

-- Input: A Grothendieck-Witt class beta over a field k
-- Output: The decomposition as a sum of hyperbolic and rank one forms

sumDecompositionString = method()
sumDecompositionString (GrothendieckWittClass) := (String) => beta -> (
    return (sumDecompositionVerbose(beta))_1
);



