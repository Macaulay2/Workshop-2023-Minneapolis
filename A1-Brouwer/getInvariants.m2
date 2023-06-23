path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");

load "diagonalize.m2"
loadPackage "RationalPoints2"
needs "GW-type.m2"
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


getHasseWittinvariant = method()

getHasseWittinvariant (GrothendieckWittClass, ZZ) := ZZ => (alpha,p) -> (
    diagAlpha := diagonalForm(alpha);
    -- Diagonalize quadratic form
    A := diagAlpha.matrix;
    n := numRows(A);
    -- Make list of diagonal entries
    -- Diagonal entries might not be integers, so we need to find an integer, n, such that when each diagonal
    -- entry turns into an integer after being multiplied by n, the lcm of the denominators of the diagonal 
    -- entries will do this, so we make a list of these denominators, and find the lcm
    diagEntries := new MutableList;
    denoDiagEntries := new MutableList;
    for i from 0 to (n-1) do(
        append(diagEntries, A_(i,i));
        append(denoDiagEntries, denominator(A_(i,i)));
    );
    -- Least common multiple of denominators
    l := lcm(denoDiagEntries);
    -- multiply the diagonal entries by the square of the lcm in order for the new diagonal entries to 
    -- represent the same quadratic form
    diagEntriesInt := (l^2)*diagEntries;

    hwi := hasseWittInvariant(diagEntriesInt,p);

    return hwi;
);


B = matrix{{1,0,2},{2,2,0},{1,0,0}}




