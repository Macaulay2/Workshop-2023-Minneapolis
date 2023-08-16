


-- Nikita Borisov and Frenly Espino

-- Given a matrix A, WittDecomp decomposes A as a sum nH + Q, where n = number of hyperbolic forms, and Q is (presumed to be) anisotropic.  

-- Input:  A square symmetric matrix A over an exact field (not RR or CC)
-- Output: (n, Q), n = integer giving number of hyperbolic spaces, Q=anisotropic forms.

-- Future Work: Should be able to input a bound.




-- WittDecomp method for InexactFieldFamily

-- WittDecompInexact calculates the Witt decomposition for matrix A over the fields kk = RR, CC
-- Function assumes that A has maximal rank 

WittDecompInexact = method()
WittDecompInexact (Matrix) := (ZZ,Matrix) => (A) -> (
    k := ring A;
    
    -- Checks that we have RR or CC as our field
    if not (instance(k,RealField) or instance(k,ComplexField) or instance(k,InexactFieldFamily)) then error "Error: base field is not RR or CC";
    
    n := numRows(A); --rank of matrix
    
    -- If k is the complex numbers, witt decomposition depends only on rank
    -- If rank is even, then matrix decomposes into n/2 hyberbolic forms with no anisotropic parts
    -- If rank is odd, matrix decomposes into (n-1)/2 hyperbolic forms with 1 x 1 anisotropic part
    if (k === CC or instance(k,ComplexField)) then (
        if (n%2 == 0) then(return (n//2,matrix(CC,{{}})))
        else return (n//2,id_(k^1));
        );
    
    -- If k is the real numbers, witt decomposition depends on rank and signature
    if (k === RR or instance(k,RealField)) then (
        diagA := congruenceDiagonalize(A);
	-- for loop counts the number of positive diagonal entries of diagA
        posEntries := 0;
	-- for loop counts the number of negative diagonal entries
        negEntries := 0;
	for i from 0 to (n - 1) do(
            if diagA_(i,i) > 0 then(
                posEntries = posEntries+1;
            );
	    if diagA_(i,i) < 0 then(
                negEntries = negEntries+1;
            );
        );

        if (posEntries + negEntries > n) then (error "A is singular");
	-- Witt index is given by how many positive-negative diagonal entry pairs exist
        wittIndex := min(posEntries,negEntries);
        signature := posEntries - negEntries; 
        if signature == 0 then (return (wittIndex,matrix(RR,{{}})))
	-- signature characterizes anisotropic part
        else if signature > 0 then ( return (wittIndex, id_(k^(signature))))
        else return (wittIndex, -id_(k^(-signature)));
        );
    );

