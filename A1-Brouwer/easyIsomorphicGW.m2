needsPackage "RationalPoints2"
needs "./GW-type.m2"

-- Tom, Andrew

------------
-- easyIsomorphicGW method
------------

-- Inputs two GW classes witnessed by matrices A and B, and tries to solve for
-- a matrix P so that P^T A P = B, exhibiting matrix congruence.
--
-- Entries of P are checked up to a certain height threshhold, hence a false
-- output from this function doesn't necessarily mean the inputs aren't isomorphic
-- in GW(k).

-- Default upper bound for height:
easyIsomorphicGW = method(
    Options => {
	HeightBound => 4  
	}
    )

easyIsomorphicGW (GrothendieckWittClass, GrothendieckWittClass) := Boolean => opts -> (beta, gamma) -> (
    -- Returns error if the matrices are defined over different base fields
    if not baseField(beta) === baseField(gamma) then error "Error: these classes have different underlying fields";
    
    A := beta.matrix;
    B := gamma.matrix;
    
    -- Return false if the matrices are not of the same size
    if not numRows(A) == numRows(B) then return false;
    
    k:= baseField(beta);
    n := numRows(A);
    
    -- Build a generic matrix P in indeterminants over our field
    R := k[x_(1,1)..x_(n,n)];
    P := genericMatrix(R,n,n);
    
    -- Take the n^2 equations produced by this matrix equality
    testZeroMatx := transpose(P)*A*P - B;
    fullEqns := flatten entries testZeroMatx;
    I := ideal(fullEqns);
    
    -- Increasing the height bound, check for solutions at each stage
    for i in 1..opts.HeightBound do(
	L := rationalPoints(k, I, Bound=>i);
	if not #L === 0 then(
	    return true;
	    break;
	    );
	if #L === 0 then(
	    print("No solutions up to height " | toString(i));
	    );
     	);
    print("Height bound exceeded");
    return false;
    )
