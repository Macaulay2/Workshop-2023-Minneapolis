--tests whether two quadratic forms are isomorphic; only implemented over RR and CC.

load "GW-type.m2"
load "diagonalize.m2"

isIsomorphic2 = method()

isIsomorphic2 (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    A=alpha.matrix;
    B=beta.matrix;
    k1=baseField(alpha);
    k2=baseField(beta);
    --base fields must be the same
    if not k1 === k2 then (
        error "Base fields are not the same";
        );
    --matrix must be square
    if numRows(A) != numColumns(A) then (
	error "Matrix is not square";
	);
    --matrix must be square
    if numRows(B) != numColumns(B) then (
	error "Matrix is not square";
	);
    --matrix must be nonsingular
    if det(A) == 0 then (
	error "Matrix is singular";
	);
    --matrix must be nonsingular
    if det(B) == 0 then (
	error "Matrix is singular";
	);
    --matrix must be symmetric
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    --matrix must be symmetric
    if B != transpose(B) then (
        error "Matrix is not symmetric";
	);
    if (k1 === CC or instance(k1,ComplexField)) then (
        if numRows(A) == numRows(B) then (
            return true;
            );
        if numRows(A) != numRows(B) then (
            return false;
           );
        )
    else if (k1 === RR or instance(k1,RealField)) then (
        diagA := diagonalize(A);
        diagB := diagonalize(B);
        posEntriesA := 0;
        posEntriesB := 0;
        for i from 0 to (numRows(A)-1) do(
            if diagA_(i,i)>0 then(
                posEntriesA=posEntriesA+1;
                );
            );
        for i from 0 to (numRows(A)-1) do(
           if diagB_(i,i)>0 then(
                posEntriesB=posEntriesB+1;
                );
            );
        if posEntriesA == posEntriesB then (
            return true;
            );
        if posEntriesA != posEntriesB then (
            return false;
            );
        )
    else print "Not implemented for fields other than RR or CC.";
    )