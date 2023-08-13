anisotropicPart = method()
anisotropicPart (Matrix) := (Matrix) => (A) -> (
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
    -- Over CC, the anisotropic part is either the rank 0 form or the rank 1 form, depending on the parity of the number of nonzero diagonal entries
    if (k === CC or instance(k,ComplexField)) then (
        nonzeroEntriesA := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) != 0 then (
                nonzeroEntriesA = nonzeroEntriesA + 1;
                );
            );
        if (anisotropicDimension(A)==0) then (
            return (diagonalMatrix(CC,{}));
            )
        else (
            return (matrix(CC,{{1}}));
            );
        )
    --Over RR, the anisotropic part consists of the positive entries in excess of the number of negative entries, or vice versa
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
        if (posEntriesA > negEntriesA) then (
            return (id_(RR^(posEntriesA-negEntriesA)));
            )
        else if (posEntriesA < negEntriesA) then (
            return (-id_(RR^(negEntriesA-posEntriesA)));
            )
        else (
            return (diagonalMatrix(RR,{}));
            );
        )
    -- Over QQ, call anisotropicPartQQ
    else if (k === QQ) then (
        error "This will call anisotropicPartQQ when that is finished"
        )
    -- Over a finite field, if the anisotropic dimension is 1, then the form is either <1> or <e>, where e is any nonsquare representative, and if the anisotropic dimension is 2 then the form is <1,-e>
    else if (instance(k, GaloisField) and k.char != 2) then (
        countNonzeroDiagA := 0;
        prodNonzeroDiagA := 1;
        for i from 0 to (numRows(A)-1) do (
	    if diagA_(i,i) != 0 then (
		countNonzeroDiagA = countNonzeroDiagA + 1;
                prodNonzeroDiagA = prodNonzeroDiagA * diagA_(i,i);
		);
	    );
        if (anisotropicDimension(A)==1) then (
            return (matrix(k,{{sub((-1)^((countNonzeroDiagA-1)/2),k)*prodNonzeroDiagA}}));
            )
        else if (anisotropicDimension(A)==0) then (
            return (diagonalMatrix(k,{}));
            )
        else (
            return (matrix(k,{{1,0},{0,sub((-1)^((countNonzeroDiagA-2)/2),k)*prodNonzeroDiagA}}));
            );
        )
    -- We should never get here
    else error "Problem with base field"
    )


anisotropicPartGW = method()
anisotropicPartGW (GrothendieckWittClass) := (GrothendieckWittClass) => (alpha) -> (
    return (gwClass(anisotropicPart(alpha.matrix)));
    )