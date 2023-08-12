---------------------------------------
-- Invariants
---------------------------------------

-- For a class defined over QQ or RR, obtain the number of positive entries in its diagonalization
numPosEntries = method()
numPosEntries (GrothendieckWittClass) := ZZ => beta ->(
    B := beta.matrix;
    n := numRows(B);
    kk := ring B;
    if not (kk === RR or instance(kk,RealField) or kk === QQ) then(
        error "Field is not QQ or RR";
        );
    diagB := congruenceDiagonalize(B);
    posEntries := 0;
    for i from 0 to (numRows(B)-1) do (
        if diagB_(i,i) > 0 then(
            posEntries = posEntries+1;
            );
	);

    return posEntries
);

-- For a class defined over QQ or RR, obtain the number of negative entries in its diagonalization
numNegEntries = method()
numNegEntries (GrothendieckWittClass) := ZZ => beta ->(
    B := beta.matrix;
    n := numRows(B);
    kk := ring B;
    if not (kk === RR or instance(kk,RealField) or kk === QQ) then(
        error "Field is not QQ or RR";
        );
    diagB := congruenceDiagonalize(B);
    negEntries := 0;
    for i from 0 to (numRows(B)-1) do (
        if diagB_(i,i) < 0 then(
            negEntries = negEntries+1;
            );
	);

    return negEntries
);

-- Returns the signature of a symmetric bilinear form over QQ or RR
signature = method()
signature (GrothendieckWittClass) := ZZ => (beta) ->(
    sig := numPosEntries(beta) - numNegEntries(beta);
    return sig
    );



---------------------------
-- Comparing forms over QQ
---------------------------

-- Inputting a form over QQ, outputs a squarefree integral representative of its discriminant
integralDiscriminant = method()
integralDiscriminant (GrothendieckWittClass) := (ZZ) => (beta) -> (
    B:= beta.matrix;
    rankForm:= numRows(B);
    kk:= ring B;
    
    if (not (kk===QQ)) then (error "GrothendieckWittClass is not over QQ");
    
    -- Take an integral diagonal representative for beta
    gamma := integralDiagonalRep(beta);
    G := gamma.matrix;
    
    discrimForm:= 1;
    for i from 0 to (rankForm-1) do(
	discrimForm = discrimForm * (G_(i,i));
	);
    
    return sub(squarefreePart(discrimForm),ZZ);
    );

-- Given a form over QQ, returns the smallest list of primes that divide its discriminat
relevantPrimes = method()
relevantPrimes (GrothendieckWittClass) := List => (beta) -> (
    B:= beta.matrix;
    rankForm:= numRows(B);
    kk:= ring B;
    
    -- Take a diagonal integral representative of the form
    gamma := integralDiagonalRep(beta);
    D := diagonalEntries(gamma);
    
    L := {};
    
    -- Append all the prime factors of each of the entries appearing on a diagonal
    for x in D do(
	L = unique(L | primeFactors(sub(x,ZZ)));
	);
    
    return L
    );



    
 -- Two Q forms over Q_p are isomorphic if they have same rank, same discriminant, and same Hasse-Witt invariant   

      
hasseWittInvariant = method()

-- epsilonHilbert computes the epsilon function for a diagonal quadratic form over Q_p
-- Function requires the list of the diagonal elements of the quadratic form, to be integers
-- Input:  A list of the diagonal elements (f_i) for the quadratic form, assumed to be integers, and a prime p
-- Output: The hasseWittInvariant function for the quadratic form (f_i) for Q_p

hasseWittInvariant (List, ZZ) := ZZ => (f,p) -> (
       a:=1;
       len:=#f;
       for i from 0 to len-1 do (
	   if not liftable(f_i,ZZ) then (error "Error:  Hilbert symbol evaluated at a non-integer");
	   );
       for i from 0 to len-2 do (
       	   for j from i+1 to len-1 do (
	       a= a * HilbertSymbol(f_i, f_j, p);
	       );
	   );
       
       return a;          
    );

hasseWittInvariant(GrothendieckWittClass, ZZ) := ZZ => (beta,p) -> (
    kk := baseField beta;
    if not (kk === QQ) then error "method is only implemented over the rationals";
    if not isPrime(p) then error "second argument must be a prime number";
    
    return hasseWittInvariant(diagonalEntries(beta),p)
    
    )




invariantFormQp =method()

-- Input: (Q, p): Quadratic form Q given by list of diagonal elements.  For now, assume list to to be integers
--                p is a prime number
-- Output:  (Rank, Disc, Hasse Invariant)

invariantFormQp (List, ZZ):= (ZZ, ZZ, ZZ) => (f, p) -> (
    -- currently will export the discriminant as a square free integer
    -- Note:  Still need a way to treat two integers as defining the same discriminant, if they differ by a
    --  square in Z_p.  They need to have the same parity of power of prime p, and the quotient of their 
    -- (prime-to-p) parts must define a unit square in Z_p^*.

    len:=#f;
    for i from 0 to (len-1) do (
	if (f_i==0) then (error "Error: Form is degenerate");
	if not liftable(f_i, ZZ ) then (error "Error: Diagonal elements of form should be integers");
	);
    a:=len;
    b:=1;
    for i from 0 to len-1 do (b=b*f_i);
    b=squarefreePart(sub(b, QQ));
    c:=hasseWittInvariant(f, p);
    return(a, b, c);
    );



-- TODO can we delete this
signatureRealQForm = method ()
signatureRealQForm (List):=(ZZ, ZZ, ZZ) => (f) -> (
    -- Input: f = list of diagonal elements of quadratic form
    -- Output: (r, s, t): Signature of form over the reals. 
    --      r= # of positive entries, s= # of negative entries, t= number of zeros
    posEntries :=0;
    negEntries:= 0;
    zeroEntries:=0;
    --for loop counts the number of negative diagonal entries
    n:=#f;
    for i from 0 to (n-1) do(
            if f_(i)>0 then(
                posEntries=posEntries+1;
            	)
	    else (
		if f_(i)<0 then(
                    negEntries=negEntries+1;
            	)
	    else (
		zeroEntries = zeroEntries +1;
		)
	    )
	);
    return (posEntries, negEntries, zeroEntries);
    
	     
    );


-- TODO can we delete?
	
-- Input: (Q): Quadratic form Q given by list of diagonal elements.  
--  	Assume list constists of integers
-- Output:  (Rank, Disc, Signature, Hasse Invariant for all primes p when not 1)	
invariantFormQ =method()
invariantFormQ (List):= (ZZ, ZZ, List, List) => (f) -> (
    -- currently will export the discriminant as a square free integer
    -- Note:  Still need a way to treat two integers as defining the same discriminant, if they differ by a
    --  square in Z_p.  They need to have the same parity of power of prime p, and the quotient of their 
    -- (prime-to-p) parts must define a unit square in Z_p^*.

    len:=#f;
    for i from 0 to (len-1) do (
	if (f_i==0) then (error "Error: Form is degenerate");
	if not liftable(f_i,ZZ ) then (error "Error: Diagonal elements of form should be integers");
	);
    a:=len;
    b:=discForm(f);
    c:=signatureRealQForm(f);
    d:=b;
    if (b<0) then (d=-b);
    -- The keys of H contain all primes dividing coefficients
    H:= hashTable( factor d);
    k:= keys H;
    l:={};
    if ((not H#?2) and hasseWittInvariant(f, 2) == -1) then l=append(l, 2);
    for i from 0 to #k-1 do (
	if  (hasseWittInvariant(f, k_i) == -1) then (
	    l=append(l,k_i);
	    );
	);
    l=sort(l);
    return (a, b, c, l);
    );

