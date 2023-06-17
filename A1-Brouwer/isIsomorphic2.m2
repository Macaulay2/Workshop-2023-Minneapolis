--Tests whether two quadratic forms are isomorphic; only implemented over QQ, RR, CC, and finite fields.

load "GW-type.m2"
load "diagonalize.m2"
load "hilbertSymbol.m2"
load "simplifyForm.m2"

isIsomorphic2 = method()

isIsomorphic2 (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    k1=baseField(alpha);
    k2=baseField(beta);
    if not (k1 === CC or instance(k1,ComplexField) or k1 === RR or instance(k1,RealField) or k1 === QQ or instance(k1, GaloisField)) then (
        error "Only implemented over QQ, RR, CC, and finite fields";
        );
    if not (k2 === CC or instance(k2,ComplexField) or k2 === RR or instance(k2,RealField) or k2 === QQ or instance(k2, GaloisField)) then (
        error "Only implemented over QQ, RR, CC, and finite fields";
        );
    A=alpha.matrix;
    B=beta.matrix;
    --Matrices must be symmetric
    if A != transpose(A) then (
        error "Underlying matrix is not symmetric";
	);
    if B != transpose(B) then (
        error "Underlying matrix is not symmetric";
	);
    if (numRows(A) != numRows(B)) then (
        return false;
        );
    diagA := diagonalize(A);
    diagB := diagonalize(B);
    if (k1 === CC or instance(k1,ComplexField)) and (k2 === CC or instance(k2,ComplexField)) then (
        nonzeroEntriesA := 0;
        nonzeroEntriesB := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) != 0 then (
                nonzeroEntriesA = nonzeroEntriesA+1;
                );
            if diagB_(i,i) != 0 then (
                nonzeroEntriesB = nonzeroEntriesB+1;
                );
            );
        if nonzeroEntriesA == nonzeroEntriesB then (
            return true;
            )
        else return false;
        )
    else if ((k1 === RR or instance(k1,RealField)) and (k2 === RR or instance(k2,RealField))) then (
        posEntriesA := 0;
        posEntriesB := 0;
        negEntriesA := 0;
        negEntriesB := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) > 0 then (
                posEntriesA=posEntriesA+1;
                );
            if diagA_(i,i) < 0 then (
                negEntriesA=negEntriesA+1;
                );
            if diagB_(i,i) > 0 then (
                posEntriesB=posEntriesB+1;
                );
            if diagB_(i,i) < 0 then (
                negEntriesB=negEntriesB+1;
                );
            );
        if ((posEntriesA == posEntriesB) and (negEntriesA == negEntriesB)) then (
            return true;
            )
        else return false;
        )
    else if ((k1 === QQ) and (k2 === QQ)) then (
        return isIsomorphicFormQ(diagA,diagB);
        )
    else if (instance(k1, GaloisField) and instance(k2, GaloisField) and k1.order == k2.order) then (
        nonzeroSquaresA := 0;
	nonzeroNonSquaresA := 0;
        nonzeroSquaresB := 0;
	nonzeroNonSquaresB := 0;
        for i from 0 to (numRows(A)-1) do (
	    if legendreBoolean(diagA_(i,i)) then (
		nonzeroSquaresA = nonzeroSquaresA + 1;
		)
	    else if diagA_(i,i) != 0 then (
		nonzeroNonSquaresA = nonzeroNonSquaresA + 1;
		);
	    if legendreBoolean(diagB_(i,i)) then (
		nonzeroSquaresB = nonzeroSquaresB + 1;
		)
	    else if diagB_(i,i) != 0 then (
		nonzeroNonSquaresB = nonzeroNonSquaresB + 1;
		);
	    );
        if ((nonzeroSquaresA == nonzeroSquaresB) and (nonzeroNonSquaresA == nonzeroNonSquaresB)) then (
            return true;
            )
        else return false;
        )
    else return false;
    )
