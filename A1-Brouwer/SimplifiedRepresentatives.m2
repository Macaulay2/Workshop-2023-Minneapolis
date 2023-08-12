------------
-- Simplifying forms
------------

---------
-- diagonalForm method
-- inputs a GWClass and outputs its diagonal form as a GWClass
---------

diagonalForm = method()
diagonalForm (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
        
    -- Check if the diagonalForm has already been computed, if so recall it from the cache
    if beta.cache.?diagonalForm then return beta.cache.diagonalForm;
    
    -- Number of rows of beta
    n := numRows(beta.matrix);
    
    -- If the field is the complex numbers, its diagonal form is an identity matrix of the same size
    if (baseField(beta) === CC or instance(baseField(beta),ComplexField)) then(
	identityMat := matrix(mutableIdentity(CC,n));
	
	-- Cache the answer before returning
	beta.cache.diagonalForm = gwClass(identityMat);
	return gwClass(identityMat)
	);
    
    -- If the field is the real numbers, we can run WittDecompInexact to determine its Witt index and anisotropic part
    if (baseField(beta) === RR or instance(baseField(beta),RealField)) then(
	
	-- Get wittIndex and anisotropic part
	(wittIndex,anisotropicPart) := WittDecompInexact(beta.matrix);
	
	-- Make empty output matrix to populate
	diagOutputMatrix := matrix(RR,{{}});
	H := matrix(RR,{{1,0},{0,-1}});
	
	
	
	-- Sum as many hyperbolic forms as the Witt index
	if wittIndex > 0 then(
	    for i in 1..wittIndex do(
		diagOutputMatrix = safeBlockSum(diagOutputMatrix,H);
		);
	    
	    );
	
	-- Add on the anisotropic part
	diagOutputMatrix = safeBlockSum(diagOutputMatrix, anisotropicPart);
	
	-- Cache and return
	beta.cache.diagonalForm = diagOutputMatrix;
	return gwClass(diagOutputMatrix)
	
	);
    
    betaMatrix := beta.matrix;
    diagonalFormOfBetaMatrix := congruenceDiagonalize(betaMatrix);
    beta.cache.diagonalForm = gwClass(diagonalFormOfBetaMatrix);
    return gwClass(diagonalFormOfBetaMatrix) 
    );


-- Extracts the diagonal entries of a matrix as a list
diagonalEntries = method()
diagonalEntries (GrothendieckWittClass) := (List) => (beta) -> (
    
    betaDiagonal := diagonalForm(beta);
    M := betaDiagonal.matrix;
    L := {};
    n := numRows M;
    
    for i from 0 to (n-1) do(
	L = append(L, M_(i,i));
	);
    return L
    );

-- Given a form over QQ, returns a diagonal form with integral entries
integralDiagonalRep = method()
integralDiagonalRep (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    kk := baseField beta;
    if not (kk === QQ) then error "method is only implemented over the rationals";
    
    L := diagonalEntries(beta);
    n := #L;
    
    -- diagonalClass takes a sequence as an input
    integralDiagonalEntries := ();
    for i from 0 to (n-1) do(
	integralDiagonalEntries = append(integralDiagonalEntries, squarefreePart(L_i))
	);
    gamma := diagonalClass(QQ,integralDiagonalEntries);
    return gamma
    );
    
    
