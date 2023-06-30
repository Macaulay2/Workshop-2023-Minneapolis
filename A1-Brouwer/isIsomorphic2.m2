-- Tests whether two quadratic forms are isomorphic; only implemented over QQ, RR, CC, and finite fields.

path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");

needs "GW-type.m2"
load "diagonalize.m2"
load "hilbertSymbol.m2"
load "simplifyForm.m2"

isIsomorphic2 = method()

isIsomorphic2 (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    k1=baseField(alpha);
    k2=baseField(beta);
    -- Ensure both base fields are supported
    if not (k1 === CC or instance(k1,ComplexField) or k1 === RR or instance(k1,RealField) or k1 === QQ or instance(k1, GaloisField)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields";
        );
    if not (k2 === CC or instance(k2,ComplexField) or k2 === RR or instance(k2,RealField) or k2 === QQ or instance(k2, GaloisField)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields";
        );
    A=alpha.matrix;
    B=beta.matrix;
    -- Ensure both underlying matrices are symmetric
    if A != transpose(A) then (
        error "Underlying matrix is not symmetric";
	);
    if B != transpose(B) then (
        error "Underlying matrix is not symmetric";
	);
    diagA := diagonalize(A);
    diagB := diagonalize(B);
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
    --Over RR, diagonal forms are equivalent if and only if they have the same number of positive, negative, and zero entries
    else if ((k1 === RR or instance(k1,RealField)) and (k2 === RR or instance(k2,RealField))) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        return ((getInvariants(alpha))_2 == (getInvariants(beta))_2);
        )
    -- Over QQ, call isIsomorphicFormQ, which checks equivalence over all completions
    else if ((k1 === QQ) and (k2 === QQ)) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        return isIsomorphicFormQ(diagA,diagB);
        )
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
