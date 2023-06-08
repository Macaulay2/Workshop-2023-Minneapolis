load "./GW-type.m2"
load "./wittDecomp.m2"

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
	(k,anisotropicPart) := wittDecompInexact(beta.matrix);
	
	-- Sum as many hyperbolic forms as the Witt index
	H := matrix(RR,{{1,0},{0,-1}});
	if k == 0 then(
	    outputMatrix := anisotropicPart;
	    );
	
	if k > 0 then(
	    outputMatrix := H;
	    for i in 1..(k-1) do(
	    	outputMatrix = outputMatrix++H;
	    	);
	    );
	if not anisotropicPart == 0 then(
	    outputMatrix = outputMatrix++anisotropicPart;
	    );
	beta.cache.diagonalForm = outputMatrix;
	return gwClass(outputMatrix)
	);
    
    A := beta.matrix;
    D := diagonalize(A);
    beta.cache.diagonalForm = gwClass(D);
    return gwClass(D) 
    )


