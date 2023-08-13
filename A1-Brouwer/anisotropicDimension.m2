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
        error "This will call anisotropicDimensionQQ when that is finished"
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
    )
