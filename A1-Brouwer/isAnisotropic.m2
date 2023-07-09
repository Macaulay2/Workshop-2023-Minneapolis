-- Tests whether a quadratic form is anisotropic; only implemented over QQ, RR, CC, and finite fields of characteristic not 2.

path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");

needs "GW-type.m2"
needs "diagonalize.m2"
needs "isAnisotropicQ.m2"

isAnisotropic = method()

isAnisotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    k=baseField(alpha);
    -- Ensure base field is supported
    if not (k === CC or instance(k,ComplexField) or k === RR or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and (k.order)%2 != 0)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    A=alpha.matrix;
    -- Ensure underlying matrix is symmetric
    if A != transpose(A) then (
        error "Underlying matrix is not symmetric";
	);
    diagA := diagonalize(A);
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
    else if (instance(k, GaloisField) and (k.order)%2 != 0) then (
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

isIsotropic = method()

isIsotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    return (not isAnisotropic(alpha));
    )
