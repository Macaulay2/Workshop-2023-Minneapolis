-- Boolean determining if a rational symmetric bilinear form is isotropic over Qp.
-- We can reduce to cases as in Serre Chapter IV Theorem 6
isIsotropicQp = method()
isIsotropicQp (GrothendieckWittClass, ZZ) := (Boolean) => (beta, p) ->(
    kk := baseField beta;
    n := numRows(beta.matrix);
    if not (kk === QQ) then error "method is only implemented over the rationals";
    if not isPrime(p) then error "second argument must be a prime number";

    -- Get the discriminant and Hasse-Witt invariant of the form at p
    d := integralDiscriminant(beta);
    
    -- Non-degenerate rank one forms are always anisotropic
    if (n==1) then(
	if det(beta.matrix) != 0 then(
	    return false;		
	    );
	);       
    
    -- If the form is rank 2, we check if d == -1 as square classes in Qp
    if (n==2) then (
	if (equalUptoPadicSquare(sub(-1, ZZ), d, p)) then (
	    return true;
	    )
	else (return false);
	);
    
    -- If the form is rank 3, we need to check if the Hilbert symbol
    -- (-1,-d)_p agrees with the Hasse-Witt invariant of the form    
    if (n==3) then(
	if (HilbertSymbol(-1, -d, p) == HasseWittInvariant(beta, p)) then (
	    return true;
	    )
	else (return false);
    	);
    
    -- If the form is rank 4, and d is a non-square, true. If d is a square and the
    -- HasseWitt invariant is equal to the Hilbert symbol (-1,-1)_p, the form is isotropic.
    -- Otherwise it isn't
    if (n==4) then(
	if ((not equalUptoPadicSquare( d, 1, p)) or (equalUptoPadicSquare(d, 1, p) and (HilbertSymbol(-1,-1,p) == HasseWittInvariant(beta, p)))) then (
	    return true;
	    )
	else (return false);		
	);    
    
    -- All forms of rank >=5 are isotropic
    if (n>=5) then(
	return true;
	);
    );

isIsotropicQp (Matrix, ZZ) := (Boolean) => (M, p) ->(
    beta := gwClass(M);
    return isIsotropicQp(beta,p)
    )    


-- We can exploit the local-to-global principle for isotropy to check 



--
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
     
  
    -- if p>2, then hilbert symbol (a,b)=1 if, a, b not divisible by p.  So HasseWittInvariant is also 1.  
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
