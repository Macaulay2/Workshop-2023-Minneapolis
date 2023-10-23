---------------------
-- Simplifying forms
---------------------

-- Input: A Grothendieck-Witt class beta
-- Output: A diagonalized form of beta

diagonalClass = method()
diagonalClass (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
        
    -- Check if the diagonalClass has already been computed, if so recall it from the cache
    if beta.cache.?diagonalClass then return beta.cache.diagonalClass;
    
    -- Number of rows of beta
    n := numRows(beta.matrix);
    
    -- If the field is the complex numbers, its diagonal form is an identity matrix of the same size
    if (baseField(beta) === CC or instance(baseField(beta),ComplexField)) then(
	identityMat := matrix(mutableIdentity(CC,n));
	
	-- Cache the answer before returning
	beta.cache.diagonalClass = gwClass(identityMat);
	return gwClass(identityMat)
	);
    
    -- If the field is the real numbers, we can diagonalize and then replace each entry on the diagon with +1 or -1 depending on its square class
    if (baseField(beta) === RR or instance(baseField(beta),RealField)) then(
	
	diagForm := congruenceDiagonalize(beta.matrix);
	L := ();
	for i from 0 to (n-1) do(
	    if diagForm_(i,i) > 0 then(
		L = append(L,1);
		);
	    if diagForm_(i,i) < 0 then(
		L = append(L,-1);
		);
	    );
	return diagonalForm(RR,L);
	
	);
    
    -- Otherwise just run congruenceDiagonalize
    betaMatrix := beta.matrix;
    diagonalClassOfBetaMatrix := congruenceDiagonalize(betaMatrix);
    
    -- The diagonal form gets cached in the GWclass type
    beta.cache.diagonalClass = gwClass(diagonalClassOfBetaMatrix);
    return gwClass(diagonalClassOfBetaMatrix) 
    );

-- Input: A GrothendieckWittClass over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: A diagonalized form of the GrothendieckWittClass, with squares stripped out

diagonalClassSimplify = method()
diagonalClassSimplify (GrothendieckWittClass) := (GrothendieckWittClass) => (alpha) -> (
    return gwClass(congruenceDiagonalizeSimplify(alpha.matrix));
    )

-- Input: A Grothendieck-Witt class beta
-- Output: The diagonal entries of beta.matrix as a list

diagonalEntries = method()
diagonalEntries (GrothendieckWittClass) := (List) => (beta) -> (
    
    betaDiagonal := diagonalClass(beta);
    M := betaDiagonal.matrix;
    L := {};
    n := numRows M;
    
    for i from 0 to (n-1) do(
	L = append(L, M_(i,i));
	);
    return L
    );

-- Input: A Grothendieck-Witt class beta over QQ
-- Output: A diagonal form with integral entries

integralDiagonalRep = method()
integralDiagonalRep (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    kk := baseField beta;
    if not (kk === QQ) then error "method is only implemented over the rationals";
    
    L := diagonalEntries(beta);
    n := #L;
    
    -- diagonalForm takes a sequence as an input
    integralDiagonalEntries := ();
    for i from 0 to (n-1) do(
	integralDiagonalEntries = append(integralDiagonalEntries, squarefreePart(L_i))
	);
    gamma := diagonalForm(QQ,integralDiagonalEntries);
    return gamma
    );
    
    
