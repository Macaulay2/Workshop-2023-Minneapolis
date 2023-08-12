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
    k1 := baseField(beta);
    k2 := baseField(gamma);
    -- Returns error if the matrices are defined over different base fields or a base field is not QQ or a finite field
    if not ((k1 === QQ  and k2 === QQ) or (instance(k1, GaloisField) and instance(k2, GaloisField) and k1.order == k2.order)) then error "These classes have non-isomorphic underlying fields, or one of the underlying fields is not QQ or a finite field";
    
    A := beta.matrix;
    B := gamma.matrix;
    
    -- Return false if the matrices are not of the same size
    if not numRows(A) == numRows(B) then return false;
    
    n := numRows(A);
    
    -- Build a generic matrix P in indeterminants over our field
    R := k1[x_(1,1)..x_(n,n)];
    P := genericMatrix(R,n,n);
    
    -- Take the n^2 equations produced by this matrix equality
    testZeroMatx := transpose(P)*A*P - substitute(B,k1);
    fullEqns := flatten entries testZeroMatx;
    I := ideal(fullEqns);
    
    -- Increasing the height bound, check for solutions at each stage
    for i in 1..(opts.HeightBound) do(
	L := rationalPoints(k1, I, Bound=>i);
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

load "isIsomorphic2.m2"
M1=matrix(GF(5),{{1,2,3},{2,4,5},{3,5,7}});
M2=matrix(GF(5),{{1,0,0},{0,3,0},{0,0,2}});
G1=gwClass(M1);
G2=gwClass(M2);
time isIsomorphic2(G1,G2);
time easyIsomorphicGW(G1,G2);

