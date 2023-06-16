path = append(path, "/home/macaulay/A1-Brouwer/");
--load "GW-type.m2"
load "wittDecomp.m2"
load "safeBlockSum.m2"
load "diagonalize.m2"
---------
-- diagonalForm method
-- inputs a GWClass and outputs its diagonal form as a GWClass
---------

diagonalForm = method()
diagonalForm (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    
    -- TODO: quick check if the form is already diagonal, if so return it and cache it
    
    
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
    
    -- If the field is the real numbers, we can run wittDecompInexact to determine its Witt index and anisotropic part
    if (baseField(beta) === RR or instance(baseField(beta),RealField)) then(
	
	-- Get wittIndex and anisotropic part
	(wittIndex,anisotropicPart) := wittDecompInexact(beta.matrix);
	
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
    diagonalFormOfBetaMatrix := diagonalize(betaMatrix);
    beta.cache.diagonalForm = gwClass(diagonalFormOfBetaMatrix);
    return gwClass(diagonalFormOfBetaMatrix) 
    );
