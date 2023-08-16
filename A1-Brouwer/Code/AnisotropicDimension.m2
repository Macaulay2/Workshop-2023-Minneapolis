---------------------------------------
-- Anisotropic dimension
---------------------------------------


-- Koprowski/Czogala Algorithm 6
isHyperbolicQp = method()
isHyperbolicQp (GrothendieckWittClass, ZZ) := Boolean => (beta, p) ->(
    B:= beta.matrix;
    rankForm:= numRows(B);
    kk:= ring B;
    
    if (not (kk===QQ)) then (error "GrothendieckWittClass is not over QQ");
    if not isPrime(p) then error "second argument must be a prime number";
    if (isDegenerate(B)) then (error "form is degenerate");
    
    -- Odd rank forms are not hyperbolic
    if odd rankForm then return false; 
    
    -- Hyperbolic forms don't have square discriminants
    d:= integralDiscriminant(beta);
    if isPadicSquare(d,p) then return false;
    
    -- At this stage, the rank and discriminant of our beta agrees with that of a hyperbolic form,
    -- so by e.g. Lam V.3.25 it suffices to check if their Hasse-Witt invariants agree
    if even rankForm then(
	
	m := sub(rankForm/2,ZZ);
	
	-- The Hasse-Witt invariant of mH:
	HasseWittHyperbolicForm := (HilbertSymbol(-1,-1,p))^(m*(m-1)/2);
	HasseWittBeta := HasseWittInvariant(beta,p);
	return (HasseWittHyperbolicForm == HasseWittBeta)
	);
    );


-- Given a symmetric bilinear form over QQ, and a prime p, outputs the anisotropic
-- dimension of the form over the p-adic numbers Q_p. Note that as all quadratic forms
-- over Q_p of dimension >=5 are isotropic, this method will always output 0, 1, 2, 3, or 4.
-- This is an implementation of Koprowski/Czogala Algorithm 8
anisotropicDimensionQp = method()
anisotropicDimensionQp (GrothendieckWittClass, ZZ) := ZZ => (beta, p) ->(
    B:= beta.matrix;
    rankForm:= numRows(B);
    kk:= ring B;
    
    if (not (kk===QQ)) then (error "GrothendieckWittClass is not over QQ");
    if not isPrime(p) then error "second argument must be a prime number";
    if (isDegenerate(B)) then (error "form is degenerate");
    
    if even rankForm then(
	-- If the form is hyperbolic it has no anisotropic part
	if isHyperbolicQp(beta,p) then return 0;
       	
	if isPadicSquare(integralDiscriminant(beta),p) then return 2;
	
	return 4;
       
	);
    
    if odd rankForm then(
	
	c := (-1)^(rankForm*(rankForm+1)/2) * integralDiscriminant(beta);
	
	gamma := gwAdd(beta, gwClass(matrix(QQ,{{c}})));
	
	if isHyperbolicQp(gamma,p) then return 1;
	
	return 3
	
	);
    );

-- Computes the anisotropic dimension of a form over QQ
-- following Algorithm 9 of Koprowski/Czogala
anisotropicDimensionQQ = method()
anisotropicDimensionQQ (GrothendieckWittClass) := ZZ => (beta) -> (
    B:= beta.matrix;
    rankForm:= numRows(B);
    kk:= ring B;
    
    if (not (kk===QQ)) then (error "GrothendieckWittClass is not over QQ");
    if (isDegenerate(B)) then (error "form is degenerate");
    
    -- The anisotropic dimension of a form over Q is the maximum of its anisotropic dimensions at any of its completions
    
    ListOfLocalAnistropicDimensions := {};
    
    -- The anisotropic dimension at RR is the absolute value of the signature of the form
    ListOfLocalAnistropicDimensions = append(ListOfLocalAnistropicDimensions, abs(signature(beta)));
    
    -- For math reasons(?) we always have to add the anisotropic dimension at the prime 2
    ListOfLocalAnistropicDimensions = append(ListOfLocalAnistropicDimensions, anisotropicDimensionQp(beta,2));
       
    -- For the remaining local fields, we can just look at relevant primes
    for p in relevantPrimes(beta) do(
	ListOfLocalAnistropicDimensions = append(ListOfLocalAnistropicDimensions, anisotropicDimensionQp(beta,p))
	
	);
    
    return max ListOfLocalAnistropicDimensions;
    );


anisotropicDimension = method()
anisotropicDimension (Matrix) := (ZZ) => (A) -> (
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
    -- Over CC, the anisotropic dimension is 0 or 1 depending on the parity of number of nonzero diagonal entries
    if (k === CC or instance(k,ComplexField)) then (
        nonzeroEntriesA := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) != 0 then (
                nonzeroEntriesA = nonzeroEntriesA + 1;
                );
            );
        return (nonzeroEntriesA%2);
        )
    --Over RR, the anisotropic dimension is the difference between the number of positive diagonal entries and the number of negative diagonal entries
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
        return (abs(posEntriesA - negEntriesA));
        )
    -- Over QQ, call anisotropicDimensionQQ
    else if (k === QQ) then (
        return anisotropicDimensionQQ(gwClass(A));
        )
    -- Over a finite field, if the number of nonzero diagonal entries is odd, then the anisotropic dimension is 1; if the number of nonzero diagonal entries is even, then the anisotropic dimension is either 0 or 2 depending on whether the nondegenerate part of the form is totally hyperbolic
    else if (instance(k, GaloisField) and k.char != 2) then (
        countNonzeroDiagA := 0;
        prodNonzeroDiagA := 1;
        for i from 0 to (numRows(A)-1) do (
	    if diagA_(i,i) != 0 then (
		countNonzeroDiagA = countNonzeroDiagA + 1;
                prodNonzeroDiagA = prodNonzeroDiagA * diagA_(i,i);
		);
	    );
        if (countNonzeroDiagA%2==1) then (
            return 1;
            )
        else if (legendreBoolean(prodNonzeroDiagA) == legendreBoolean(sub((-1)^(countNonzeroDiagA/2),k))) then (
            return 0;
            )
        else (
            return 2;
            );
        )
    -- We should never get here
    else error "Problem with base field"
    )


anisotropicDimension (GrothendieckWittClass) := (ZZ) => (alpha) -> (
    return(anisotropicDimension(alpha.matrix));
    );

isotropicDimension = method()
isotropicDimension (GrothendieckWittClass) := ZZ -> (alpha) -> (
    n := numRows(alpha.matrix);
    return (n - anisotropicDimension(alpha))
    );

WittIndex = method()
WittIndex (GrothendieckWittClass) := ZZ -> (alpha) -> (
    return isotropicDimension(alpha)
    );
