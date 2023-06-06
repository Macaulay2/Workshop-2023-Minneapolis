load "./GW-type.m2"
load "./wittDecomp.m2"
diagonalForm = method()

diagonalForm (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    
    n:= numRows(beta.matrix);
    if (baseField(beta) === CC or instance(baseField(beta),ComplexField)) then(
	identityMat := matrix(mutableIdentity(CC,n));
	return gwClass(identityMat)
	);
    
    if (baseField(beta) === RR or instance(baseField(beta),RealField)) then(
	(k,anisotropicPart) := wittDecompInexact(beta.matrix);
	D:= matrix(RR,{{1,0},{0,-1}});
	if k>1 then(
	    for i in 1..(k-1) do(
	    	D = D++D;
	    	);
	    );
	outputMatrix:=D;
	if not anisotropicPart == 0 then(
	    outputMatrix = outputMatrix++anisotropicPart;
	    );
	return outputMatrix
	);
    
    A := beta.matrix;
    D := diagonalize(A);
    return gwClass(D)
    )
