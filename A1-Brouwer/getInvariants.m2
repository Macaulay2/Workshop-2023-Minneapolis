load "./diagonalize.m2"
loadPackage "RationalPoints2"
load "GW-type.m2"
-- getInvariants method spits out basic invariants of a matrix representing a bilinear form
-- isIdenticalDiscriminant tells us whether two matrices representing bilinear forms have the same discriminant

-- Authors: Frenly, Andrew

getInvariants = method()
getInvariants (GrothendieckWittClass) := List => (alpha) ->(
    A := alpha.matrix;
    n := numRows(A);
    kk := ring A;

    --Representative of discriminant of bilinear form
    discrRep := det A;

    if (kk === RR or instance(kk,RealField) or kk === QQ) then(
        diagA := diagonalize(A);
        posEntries := 0;
        negEntries := 0;
	--For-loop counts the number of positive, negative, and zero diagonal entries of diagA
        zeroEntries := 0;
         for i from 0 to (n-1) do(
            if diagA_(i,i) > 0 then(
                posEntries = posEntries+1;
            );
            if diagA_(i,i) == 0 then(
                zeroEntries = zeroEntries+1;
            );
            if diagA_(i,i) < 0 then(
                negEntries = negEntries+1;
            );
        );

    	-- Signature is a list
        signature := {posEntries,zeroEntries,negEntries};

    	-- Signature is the 2nd entry of output list
        return {n,signature,discrRep});

    -- Do not output signature if not working over ordered fields
    {n,discrRep}
    )

isIdenticalDiscriminant = method()

isIdenticalDiscriminant (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    A := alpha.matrix;
    B := beta.matrix;
    a := det(A);
    b := det(B);
    kk := ring A;
    R := kk[x];

    -- This will have roots in field kk if and only if b and a represent the same discriminant
    f = (x^2) - (a/b);

    L := zeros f;

    -- Return false if roots don't exist in kk
    if (L == {}) then(return false)
    else return true;
);

B = matrix{{1,0,2},{2,2,0},{1,0,0}}




