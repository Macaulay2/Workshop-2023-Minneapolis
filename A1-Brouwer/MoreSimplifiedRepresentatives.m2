-- Diagonalizes a Grothendieck--Witt class
diagonalFormNoSimplify = method()
diagonalFormNoSimplify (GrothendieckWittClass) := (GrothendieckWittClass) => (alpha) -> (
    return gwClass(congruenceDiagonalize(alpha.matrix));
    )

-- Diagonalizes a Grothendieck--Witt class over QQ, RR, CC, or a finite field of characteristic not 2 and strips out squares
diagonalFormSimplify = method()
diagonalFormSimplify (GrothendieckWittClass) := (GrothendieckWittClass) => (alpha) -> (
    return gwClass(congruenceDiagonalizeSimplify(alpha.matrix));
    )

-- Diagonalizes a quadratic form over QQ, RR, CC, or a finite field of characteristic not 2 and strips out squares
congruenceDiagonalizeSimplify = method()
congruenceDiagonalizeSimplify (Matrix) := (Matrix) => (AnonMut) -> (
    k := ring AnonMut;
    if not (k === CC or instance(k,ComplexField) or k === RR or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and k.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    if not isSquareAndSymmetric(AnonMut) then error "matrix is not symmetric";

    diagForm := congruenceDiagonalize(AnonMut);
    A := mutableMatrix diagForm;
    n := numRows(A);

    -- If the field is the complex numbers, we can replace each nonzero entry of the diagonalization by 1
    if (k === CC or instance(k,ComplexField)) then (
	for i from 0 to (n-1) do (
	    if A_(i,i) != 0 then (
		A_(i,i) = 1;
		);
            );
	return (matrix A);
	)
    
    -- If the field is the real numbers, we can replace each positive entry of the diagonalization by 1 and each negative entry by -1
    else if (k === RR or instance(k,RealField)) then (
	for i from 0 to (n-1) do (
	    if A_(i,i) > 0 then (
		A_(i,i) = 1;
		);
	    if A_(i,i) < 0 then (
		A_(i,i) = -1;
		);
	    );
	return (matrix A);
	)

    -- If the field is the rational numbers, we can diagonalize and take squarefree parts
    else if (k === QQ) then (
	for i from 0 to (n-1) do (
            A_(i,i) = squarefreePart(A_(i,i));
	    );
	return (matrix A);
	)

    -- Over a finite field, we can diagonalize and replace every entry by 1 or a nonsquare representative
    else if (instance(k, GaloisField) and k.char != 2) then (
        nonSquareRep := sub(-1,k);
        if (legendreBoolean(sub(-1,k))) then (
	    for i from 0 to (n-1) do (
	        if (diagForm_(i,i) != 0 and (not legendreBoolean(diagForm_(i,i)))) then (
	     	    nonSquareRep = diagForm_(i,i);
                    break;
		    );
                );
	    );
	for i from 0 to (n-1) do (
	    if (A_(i,i) != 0 and legendreBoolean(A_(i,i))) then (
		A_(i,i) = 1;
		);
	    if (A_(i,i) != 0 and not legendreBoolean(A_(i,i))) then (
		A_(i,i) = nonSquareRep;
		);
	    );
	return (matrix A);
	)

    -- We should never get here
    else error "Problem with base field"
    )
