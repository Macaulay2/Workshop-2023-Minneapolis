-- Boolean checking if two forms over Q are isomorphic
isIsomorphicFormQ = method()
isIsomorphicFormQ (GrothendieckWittClass, GrothendieckWittClass) := Boolean => (alpha, beta) -> (
    if (not baseField(alpha) === QQ) then error "first input must be a form defined over QQ";
    if (not baseField(beta) === QQ) then error "second input must be a form defined over QQ";
    
    -- Check if the ranks agree
    if (not numRows(alpha.matrix) == numRows(beta.matrix)) then return false;
    
    -- Check if the signatures (Hasse-Witt invariants at RR) agree
    if (not signature(alpha) == signature(beta)) then return false;
    
    -- Check if the discriminants agree
    if (not integralDiscriminant(alpha) == integralDiscriminant(beta)) then return false;
    
    -- Check if all the Hasse-Witt invariants agree
    PrimesToCheck := unique(relevantPrimes(alpha) | relevantPrimes(beta));
    flag := 0;
    for p in PrimesToCheck do(
	if (HasseWittInvariant(alpha,p) != HasseWittInvariant(beta,p)) then(
	    flag = 1;
	    break;
	    );
	);
    if flag == 0 then (
	return true;
	)
    else (
	return false;
	);
    );

isIsomorphicFormQ (Matrix, Matrix) := Boolean => (M,N) -> (
    return isIsomorphicFormQ(gwClass(M),gwClass(N));
    );



-- Boolean checking if two Grothendieck-Witt classes are the same
gwIsomorphic = method()
gwIsomorphic (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    k1:=baseField(alpha);
    k2:=baseField(beta);
    -- Ensure both base fields are supported
    if not (k1 === CC or instance(k1,ComplexField) or k1 === RR or instance(k1,RealField) or k1 === QQ or (instance(k1, GaloisField) and k1.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    if not (k2 === CC or instance(k2,ComplexField) or k2 === RR or instance(k2,RealField) or k2 === QQ or (instance(k1, GaloisField) and k1.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    
    A:=alpha.matrix;
    B:=beta.matrix;
    
    -- Ensure both underlying matrices are symmetric
    if not isSquareAndSymmetric(A) then error "Underlying matrix is not symmetric";
    if not isSquareAndSymmetric(B) then error "Underlying matrix is not symmetric";
    
    diagA := congruenceDiagonalize(A);
    diagB := congruenceDiagonalize(B);
    
    -----------------------------------
    -- Complex numbers
    -----------------------------------
    
    -- Over CC, diagonal forms over spaces of the same dimension are equivalent if and only if they have the same number of nonzero entries
    if (k1 === CC or instance(k1,ComplexField)) and (k2 === CC or instance(k2,ComplexField)) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        nonzeroEntriesA := 0;
        nonzeroEntriesB := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) != 0 then (
                nonzeroEntriesA = nonzeroEntriesA + 1;
                );
            if diagB_(i,i) != 0 then (
                nonzeroEntriesB = nonzeroEntriesB + 1;
                );
            );
        return (nonzeroEntriesA == nonzeroEntriesB);
        )
    
    -----------------------------------
    -- Real numbers
    -----------------------------------
    
    --Over RR, diagonal forms are equivalent if and only if they have the same number of positive, negative, and zero entries
    else if ((k1 === RR or instance(k1,RealField)) and (k2 === RR or instance(k2,RealField))) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        return (signature(alpha)==signature(beta));
        )
    
    -----------------------------------
    -- Rational numbers
    -----------------------------------
    
    -- Over QQ, call isIsomorphicFormQ, which checks equivalence over all completions
    else if ((k1 === QQ) and (k2 === QQ)) then (
        return isIsomorphicFormQ(alpha,beta);
        )
    
    -----------------------------------
    -- Finite fields
    -----------------------------------
    
    -- Over a finite field, diagonal forms over spaces of the same dimension are equivalent if and only if they have the same number of nonzero entries and the product of these nonzero entries is in the same square class
    else if (instance(k1, GaloisField) and instance(k2, GaloisField) and k1.order == k2.order) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        countNonzeroDiagA := 0;
        countNonzeroDiagB := 0;
        prodNonzeroDiagA := 1;
        prodNonzeroDiagB := 1;
        for i from 0 to (numRows(A)-1) do (
	    if diagA_(i,i) != 0 then (
		countNonzeroDiagA = countNonzeroDiagA + 1;
                prodNonzeroDiagA = prodNonzeroDiagA * diagA_(i,i);
		);
	    if diagB_(i,i) != 0 then (
		countNonzeroDiagB = countNonzeroDiagB + 1;
                prodNonzeroDiagB = prodNonzeroDiagB * diagB_(i,i);
		);
	    );
        return ((countNonzeroDiagA == countNonzeroDiagB) and (legendreBoolean(prodNonzeroDiagA) == legendreBoolean(prodNonzeroDiagB)));
        )
    -- If we get here, the base fields are not isomorphic
    else error "Base fields are not isomorphic"
    )

-- Boolean checking whether two symmetric bilinear forms are isometric
isIsometricForm = method()
isIsometricForm (Matrix,Matrix) := (Boolean) => (A,B) -> (
    k1:=ring A;
    k2:=ring B;
    -- Ensure both base fields are supported
    if not (k1 === CC or instance(k1,ComplexField) or k1 === RR or instance(k1,RealField) or k1 === QQ or (instance(k1, GaloisField) and k1.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    if not (k2 === CC or instance(k2,ComplexField) or k2 === RR or instance(k2,RealField) or k2 === QQ or (instance(k1, GaloisField) and k1.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    
    -- Ensure both matrices are symmetric
    if not isSquareAndSymmetric(A) then error "Underlying matrix is not symmetric";
    if not isSquareAndSymmetric(B) then error "Underlying matrix is not symmetric";
    
    diagA := congruenceDiagonalize(A);
    diagB := congruenceDiagonalize(B);
    
    -----------------------------------
    -- Complex numbers
    -----------------------------------
    
    -- Over CC, diagonal forms over spaces of the same dimension are equivalent if and only if they have the same number of nonzero entries
    if (k1 === CC or instance(k1,ComplexField)) and (k2 === CC or instance(k2,ComplexField)) then (
        return ((numRows(A) == numRows(B)) and (numNonzeroDiagEntries(diagA) == numNonzeroDiagEntries(diagB)));
        )
    
    -----------------------------------
    -- Real numbers
    -----------------------------------
    
    --Over RR, diagonal forms of the same dimension are equivalent if and only if they have the same number of positive and negative entries
    else if ((k1 === RR or instance(k1,RealField)) and (k2 === RR or instance(k2,RealField))) then (
        return ((numRows(A) == numRows(B)) and (numPosDiagEntries(diagA) == numPosDiagEntries(diagB)) and (numNegDiagEntries(diagA) == numNegDiagEntries(diagB)));
        )
    
    -----------------------------------
    -- Rational numbers
    -----------------------------------
    
    -- Over QQ, if spaces have same dimension and nondegenerate parts have same dimension, then call isIsomorphicFormQ, which checks equivalence over all completions
    else if ((k1 === QQ) and (k2 === QQ)) then (
        return ((numRows(A) == numRows(B)) and (numRows(nondegeneratePartDiagonal(A)) == numRows(nondegeneratePartDiagonal(B))) and (isIsomorphicFormQ(gwClass(nondegeneratePartDiagonal(A)),gwClass(nondegeneratePartDiagonal(B)))));
        )
    
    -----------------------------------
    -- Finite fields
    -----------------------------------
    
    -- Over a finite field, diagonal forms over spaces of the same dimension are equivalent if and only if they have the same number of nonzero entries and the product of these nonzero entries is in the same square class
    else if (instance(k1, GaloisField) and instance(k2, GaloisField) and k1.char !=2 and k2.char != 2 and k1.order == k2.order) then (
        return (numRows(A) == numRows(B)) and (numNonzeroDiagEntries(diagA) == numNonzeroDiagEntries(diagB)) and (det(nondegeneratePartDiagonal(A)) == sub(det(nondegeneratePartDiagonal(B)),k1));
        )
    -- If we get here, the base fields are not isomorphic
    else error "Base fields are not isomorphic"
    )
