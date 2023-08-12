

-- Boolean checking if a diagonal form over QQ with integer coefficients is isotropic locally at a prime p
isIsotropicDiagFormQp = method()
isIsotropicDiagFormQp (List, ZZ):=(Boolean) => (f, p) -> (
    
    n:=#f;
    d:=lift(discForm(f),ZZ);
    e:=hasseWittInvariant(f, p);
    if (n==2) then (
    -- need to compare d ==-1 in Qp^*/Qp^2	
	if (equalUptoPadicSquare(sub(-1, ZZ), d, p)) then (
	    return true;
	    )
	else (return false);
	)
    else (
	if (n==3) then (
	    if (HilbertSymbol(-1, -d, p) == hasseWittInvariant(f, p)) then (
		return true;
		)
	    else (return false);
	    )
	else (
	    if (n==4) then (
		if ((not equalUptoPadicSquare( d, 1, p)) or (equalUptoPadicSquare(d, 1, p) and (HilbertSymbol(-1,-1,p) == hasseWittInvariant(f, p)))) then (
		    return true;
		    )
		else (return false);
		)
	    else (
		return true;
		)
	    )
	)
    );




-- Boolean determining if a form is isotropic over Qp
isIsotropicQp = method()
isIsotropicQp (GrothendieckWittClass, ZZ) := (Boolean) => (beta, p) ->(
    kk := baseField beta;
    if not (kk === QQ) then error "method is only implemented over the rationals";
    if not isPrime(p) then error "second argument must be a prime number";
    
    -- Take a diagonal integral representative of beta
    betaInt := integralDiagonalRep(beta);
    
    -- Extract its diagonal entries as a list
    L := toList(diagonalEntries(betaInt));
    
    -- Call the previous function
    return isIsotropicDiagFormQp(L,p)
    )

isIsotropicQp (Matrix, ZZ) := (Boolean) => (M, p) ->(
    beta := gwClass(M);
    return isIsotropicQp(beta,p)
    )    
    


-- Boolean checking if a diagonal form over QQ with integer coefficients is isotropic at every local field including RR
isIsotropicDiagFormQ = method()
isIsotropicDiagFormQ (List):=(Boolean) => (f) ->(
    -- need to check if anisotropic over RR and all Qp
    n:= #f;
    a:=signatureRealQForm(f);
    -- if real form is definite then form f is not isotropic
    if (a_0 * a_1 <=0 ) then (
	return false;
	) 
    else (
	if (n>4) then (return true; ) 
	else (
	    if (n==1 ) then (return false; );
	    -- for n=2,3,4, the form is isotropic if p!=2 and p doesn't divide any of the 
	    -- diagonal terms as the conditions in Q_p are automatically satisfied.
	    d:= discForm(f);
	    d1:= d;
	    if (d1<0) then (d1=-d1);
	    H:= hashTable( factor d1);
	    k:= keys H;
	    i:=0;
	    l:= #k;
	    if (not isIsotropicDiagFormQp(f, 2)) then (return false);
	    while (i< l) do (
	    	p := k_i;
		if ( not isIsotropicDiagFormQp(f, p)) then (return false);
		i=i+1;
		);
	    return true;
	    );
	);
    );

	


--------------------------
-- Checking isotropy
--------------------------

-- isAnisotropicDiagFormQp determines if a diagonal quadratic form with integral coefficients is anisotropic over the p-adic field Q_p
-- Methods cited here can be founded in Serre's "A course in arithmetic," Chapter III, Theorem 6
  
isAnisotropicDiagQp = method()
isAnisotropicDiagQp (List, ZZ) := (Boolean) => (f, p) -> (
    -- Input: (f, p): f list of integrers (the diagonal elements of the form), p an integral prime
    -- Output: true if form does not represent 0 over Qp
    
    -- Check if f is list, p is prime 
    if not isPrime(p) then error "Error: isAnisotropicDiagFormQp called with p not an integer prime";
    n:=#f;
    for i from 0 to n-1 do (
    	if (not  liftable(f_i,ZZ)) then (error "Error: isAnisotropicDiagFormQp called with a non-integral list");
    );
    d:=discForm(f);
    
    -- if discriminant is 0, then form is degenerate and is isotropic
    if (d==0) then (return false);
    
    -- can now assume form is nondegenerate
    -- a rank 1 non-degenerate form is anisotropic
    if (n==1) then (return true);
    -- a rank >=5 form is isotropic
    if (n>4) then (return false);
    
    -- now only need to consider ranks 2, 3, 4
    
    e:=hasseWittInvariant(f, p);
    
    if (n==2) then (
    -- need to compare d ==-1 in Qp^*/Qp^2; if so, form is isotropic	
	if (equalUptoPadicSquare(sub(-1, ZZ), d, p)) then (
	    return false;
	    )
	else (return true);
	)
    else (
	-- now use the criteria for n=3
	-- need to check if (-1,-d)=hasseWittInvariant for f ; if equal, form is isotropic
	if (n==3) then (
	    if (HilbertSymbol(-1, -d, p) == hasseWittInvariant(f, p)) then (
		return false;
		)
	    else (return true);
	    )
	else (
	 -- now use the criteria for n=4
	-- need to check (d=1 and (-1,-1) != =hasseWittInvariant for f) for form to be anisotropic 
	    if (n==4) then (
		if ((equalUptoPadicSquare( d, 1, p)) and  (not (HilbertSymbol(-1,-1,p) == hasseWittInvariant(f, p)))) then (
		    return true;
		    )
		else (return false);
		)
	    
	    -- TODO - we should be able to delete this else case since it won't ever be hit?
	    else (
		error "Error: rank should have been 2, 3, 4, but isn't";
		return false;
		)
	    )
	)
    );


--isAnisotropicQ takes n GWClass over QQ and returns a Boolean based on whether or not 
--the class is anisotropic
--unlike sumDecomposition it can say for certain if the form is anisotropic or not since it
--does not use rationalPoints

-- Input:  A GrothendieckWittClass for a quadratic form
-- Output: True if form is Anisotropic;  False, if form is Isotropic


isAnisotropicQ = method()
isAnisotropicQ (GrothendieckWittClass) := Boolean => (alpha) -> (
    A:= alpha.matrix;
    n:= numRows(A);
    kk:= ring A;

    if (not (kk===QQ)) then (error "GrothendieckWittClass is not over QQ");

    --check if form is degenerate
    if (isDegenerate(A)) then (return false);
    
     -- if rank =1, then a non-degenerate form is anisotropic
    if (n==1) then (return true);
    
    --if rank>=5, we can use signature do decide this
    --the non-degenerate form will be anisotropic iff all diagonal entries have same sign
    if (n>= 5) then (
        return ((numPosEntries(alpha) == n) or (numNegEntries(alpha) == n));
    );
   

    --if 2<= rank <=4, we need to take p-adic completions
    -- First, we diagonalize the matrix
    diagA := congruenceDiagonalize(A);  
    -- Then obtain the diagonal entries.  These will be rational numbers;
    diagEntriesA := apply(n, i-> A_(i,i));
    -- Make then integers by multiplying by integer squares (so that the forms are equivalent);
    -- The sub command forces the list to be integers
    diagIntEntriesA:= apply(n, i-> squarefreePart(sub(numerator(diagEntriesA_i) * denominator(diagEntriesA_i),ZZ)));
    
    -- disc = discriminant of form, product of diagonal elements.  
    disc:= discForm(diagIntEntriesA);
    
    
   
   
    -- Using Q_p criteria from Thm 6, Section 2.2 of Serre's Course in Arithmetic
   
    -- For n=2, need -disc to be a square for all p to be isotropic, so in particular, need -disc=1 for isotropic
    if (n==2) then (  
	 -- Make disc a squarefree integer
	 d2 := squarefreePart(disc); 
	 return (not (d2==-1) ) 
	 );
     
  
    -- if p>2, then hilbert symbol (a,b)=1 if, a, b not divisible by p.  So hasseWittInvariant is also 1.  
    -- Then for n=3,4, form is automatically isotropic if p doesn't divide disc. 
    -- So only need to check if f is anisotropic over Q_p for p=2 and primes p dividing disc.
    
    -- first check p=2 case
    if (isAnisotropicDiagQp(diagIntEntriesA, 2)) then (return true);
    
    -- create list of primes dividing disc
    -- first take absolute value of disc
    d1:= disc;
    if (d1<0) then (d1=-d1);
    
    -- H is HashTable of factors of disc
    
    H:= hashTable( factor d1);
    -- the keys k are the prime factors
    k:= keys H;
    i:=0;
   
    while (i< #k) do (
	    	p := k_i;
		if (isAnisotropicDiagQp(diagIntEntriesA, p)) then (return true);
		i=i+1;
		);

-- if the function hasn't returned false yet, then isotropic over all primes p.  hence form is isotropic over Q	    
  
    return false;
);
	

-- Boolean returning if a symmetric bilinear form is anisotropic	
isAnisotropic = method()
isAnisotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    k:=baseField(alpha);
    -- Ensure base field is supported
    if not (k === CC or instance(k,ComplexField) or k === RR or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and k.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    A:=alpha.matrix;
    -- Ensure underlying matrix is symmetric
    if A != transpose(A) then (
        error "Underlying matrix is not symmetric";
	);
    diagA := congruenceDiagonalize(A);
    -- Over CC, a diagonal form is anisotropic if and only if it is nondegenerate and has dimension 0 or 1
    if (k === CC or instance(k,ComplexField)) then (
        nonzeroEntriesA := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) != 0 then (
                nonzeroEntriesA = nonzeroEntriesA + 1;
                );
            );
        return (nonzeroEntriesA == numRows(A) and numRows(A) <= 1);
        )
    --Over RR, a diagonal form is anisotropic if and only if all of its diagonal entries are positive or all of its diagonal entries are negative
    else if (k === RR or instance(k,RealField)) then (
        posEntriesA := 0;
        negEntriesA := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) > 0 then (
                posEntriesA = posEntriesA + 1;
                );
            if diagA_(i,i) < 0 then (
                negEntriesA = negEntriesA + 1;
                );
            );
        return ((posEntriesA == numRows(A)) or (negEntriesA == numRows(A)));
        )
    -- Over QQ, call isAnisotropicQ
    else if (k === QQ) then (
        return isAnisotropicQ(alpha);
        )
    -- Over a finite field, a diagonal form is anisotropic if and only if it is nondegenerate, of dimension at most 2, and not the hyperbolic form 
    else if (instance(k, GaloisField) and k.char != 2) then (
        countNonzeroDiagA := 0;
        prodNonzeroDiagA := 1;
        for i from 0 to (numRows(A)-1) do (
	    if diagA_(i,i) != 0 then (
		countNonzeroDiagA = countNonzeroDiagA + 1;
                prodNonzeroDiagA = prodNonzeroDiagA * diagA_(i,i);
		);
	    );
        return ((countNonzeroDiagA == numRows(A)) and (numRows(A) <= 1 or (numRows(A) == 2 and  legendreBoolean(prodNonzeroDiagA) != legendreBoolean(sub(-1,k)))));
        )
    -- We should never get here
    else error "Problem with base field"
    )

-- Boolean returning if a symmetric bilinear form is isotropic	
isIsotropic = method()
isIsotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    return (not isAnisotropic(alpha));
    )
