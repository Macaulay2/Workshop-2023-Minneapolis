isIsomorphicDiagFormQ = method ()
isIsomorphicDiagFormQ (List, List):= (Boolean) => (f, g) -> (
    -- Input: (f,g):  f, g are lists of diagonal elements (in integers) of two forms over Q
    --         
    -- Output: true if f, g are isomorphic over Q
    
    -- Method:  Calculate rank, discriminant over Q, signature over R
    --          These must all agree.
    --	      	If so, then need to see if Hasse invariant agrees for all primes p.
    --          Hasse symbol is 1 if prime p doesn't divide any of the diagonal elts of form
    --          So (i) determine the list of primes dividing any of the diagonal elts. 
    --             (ii) Calculate and compare Hasse invariant at each p.  If equal for all p,
    --                  then isomorphic forms

    -- Compare signatures over R
    disc1 := squarefreePart(discForm(f));
    disc2 := squarefreePart(discForm(g));
    
    if (signatureRealQForm(f) == signatureRealQForm(g) and disc1 == disc2) then (
	    -- signature and discriminants agree. Now need to test Hasse invariant
	    d:= disc1 * disc2;
	    if (d<0) then (d = -d);
	    if (d==0) then error "Error: Form is degenerate";
	    H:=hashTable (factor d);
	    k:= keys H;    
	    i:=0;
	    flag:=0;
	    while (i< #k and flag==0)  do (
		p:=k_i;
		if (hasseWittInvariant(f, p) != hasseWittInvariant(g, p)) then (
		    flag=1;
		    );
		i=i+1;
		);
    	    if flag==0 then (
		return true;
		)
    	    else (
		return false;
		);
	    )
	else (
	    return false;
	    );
    );




-- Boolean checking if two symmetric bilinear forms over QQ are isomorphic
isIsomorphicFormQ = method ()
isIsomorphicFormQ (Matrix, Matrix):= (Boolean) => (f, g) -> (
    
    -- Check same size
    if (numRows f != numRows g) then (return false;);
    
    -- First, we diagonalize both matrices
    df:= congruenceDiagonalize f;
    dg:= congruenceDiagonalize g;
    
   
    -- Then make all entries integers, by multiplying by denominator squared, and clearing squares
    -- create list of diagonal entries in integers
    n:=numRows f;
    f1:=apply(n, i-> df_(i,i));
    g1:=apply(n, i-> dg_(i,i));
    
  
    -- convert rational diagonals to square-free integers by multiplying by squares
    
    f2:= apply(n, i-> squarefreePart(sub(numerator(f1_i) * denominator(f1_i),ZZ)));
    g2:= apply(n, i-> squarefreePart(sub(numerator(g1_i) * denominator(g1_i),ZZ)));
    
    -- Now compare forms
    return isIsomorphicDiagFormQ(f2, g2);
    
    );
    




-- Boolean checking if two Grothendieck-Witt classes are the same
gwIsomorphic = method()
gwIsomorphic (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    k1:=baseField(alpha);
    k2:=baseField(beta);
    -- Ensure both base fields are supported
    if not (k1 === CC or instance(k1,ComplexField) or k1 === RR or instance(k1,RealField) or k1 === QQ or instance(k1, GaloisField)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields";
        );
    if not (k2 === CC or instance(k2,ComplexField) or k2 === RR or instance(k2,RealField) or k2 === QQ or instance(k2, GaloisField)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields";
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
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        return isIsomorphicFormQ(diagA,diagB);
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
