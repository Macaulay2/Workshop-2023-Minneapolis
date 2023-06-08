--tests whether two quadratic forms are isomorphic; only implemented over RR and CC.

load "GW-type.m2"
load "diagonalize.m2"

isIsomorphic2 = method()

isIsomorphic2 (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    A=alpha.matrix;
    B=beta.matrix;
    k1=baseField(alpha);
    k2=baseField(beta);
    --matrix must be symmetric
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    --matrix must be symmetric
    if B != transpose(B) then (
        error "Matrix is not symmetric";
	);
    if (k1 === CC or instance(k1,ComplexField)) and (k2 === CC or instance(k2,ComplexField)) then (
        diagA := diagonalize(A);
        diagB := diagonalize(B);
        nonzeroEntriesA := 0;
        nonzeroEntriesB := 0;
        for i from 0 to (numRows(A)-1) do(
            if diagA_(i,i) != 0 then(
                nonzeroEntriesA = nonzeroEntriesA+1;
                );
            );
        for i from 0 to (numRows(A)-1) do(
           if diagB_(i,i) != 0 then(
                nonzeroEntriesB = nonzeroEntriesB+1;
                );
            );
        if nonzeroEntriesA == nonzeroEntriesB then (
            return true;
            )
        else return false;
        )
    else if ((k1 === RR or instance(k1,RealField)) and (k2 === RR or instance(k2,RealField))) then (
        diagA := diagonalize(A);
        diagB := diagonalize(B);
        posEntriesA := 0;
        posEntriesB := 0;
        negEntriesA := 0;
        negEntriesB := 0;
        for i from 0 to (numRows(A)-1) do(
            if diagA_(i,i)>0 then(
                posEntriesA=posEntriesA+1;
                );
            if diagA_(i,i)<0 then(
                negEntriesA=negEntriesA+1;
                );
            );
        for i from 0 to (numRows(A)-1) do(
           if diagB_(i,i)>0 then(
                posEntriesB=posEntriesB+1;
                );
           if diagB_(i,i)<0 then(
                negEntriesB=negEntriesB+1;
                );
            );
        if ((posEntriesA == posEntriesB) and (negEntriesA == negEntriesB)) then (
            return true;
            )
        else return false;
        )
    else error "Only implemented if both base fields are RR or both base fields are CC";
    )
