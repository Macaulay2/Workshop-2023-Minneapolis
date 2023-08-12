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
    if isPadicSquare(d) then return false;
    
    -- At this stage, the rank and discriminant of our beta agrees with that of a hyperbolic form,
    -- so by e.g. Lam V.3.25 it suffices to check if their Hasse-Witt invariants agree
    if even rankForm then(
	
	m := sub(rankForm/2,ZZ);
	
	-- The Hasse-Witt invariant of mH:
	HasseWittHyperbolicForm := HilbertSymbol(-1,-1,p)^{m*(m-1)/2};
	HasseWittBeta := hasseWittInvariant(beta,p);
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
       	
	if isPadicSquare(integralDiscriminant(beta)) then return 2;
	
	return 4;
       
	);
    
    if odd rankForm then(
	
	c := (-1)^{rankForm*(rankForm+1)/2} * integralDiscriminant(beta);
	
	gamma := gwAdd(beta, gwClass(QQ,{{c}}));
	
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


-- -- Computes the anistropic dimension of a form over a finite field
-- anisotropicDimensionGF = method()
-- anisotropicDimension (GrothendieckWittClass) := ZZ => (beta) -> (
    
    
--     )
