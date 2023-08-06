-- Tom

path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");

needs "discForm.m2"
needs "HilbertSymbol.m2"


-- isAnisotropicDiagFormQp determines if a diagonal quadratic form with integral coefficients is anisotropic
  
isAnisotropicDiagQp = method()

isAnisotropicDiagQp (List, ZZ) := (Boolean) => (f, p) -> (
    -- Input: (f, p): f list of integrers (the diagonal elements of the form), p an integral prime
    -- Output: true if form does not represent 0 over Qp
    
    -- Check if f is list, p is prime 
    if (not isPrime(p) or (not (ring p === ZZ))) then (print "Error: isAnisotropicDiagFormQp called with p not an integer prime");
    n:=#f;
    for i from 0 to n-1 do (
    	if (not (ring f_i === ZZ)) then (print "Error: isAnisotropicDiagFormQp called with a non-integral list");
    );
    d:=discForm(f);
    
    -- if discriminant is 0, then form is degenerate and is isotropic
    if (d==0) then (return false);
    
    -- can now assume form is nondegenerate
    -- a rank 1 non-degenerate form is anisotropic
    if (n==1) then (return true);
    -- a rank >=5 form is isotropic
    if (n>4) then (return false);
    
    -- now only need to consider ranks 2, 3, 4
    
    e:=hasseWittInvariant(f, p);
    
    if (n==2) then (
    -- need to compare d ==-1 in Qp^*/Qp^2; if so, form is isotropic	
	if (equalUptoPadicSquare(sub(-1, ZZ), d, p)) then (
	    return false;
	    )
	else (return true);
	)
    else (
	-- now use the criteria for n=3
	-- need to check if (-1,-d)=hasseWittInvariant for f ; if equal, form is isotropic
	if (n==3) then (
	    if (HilbertSymbol(-1, -d, p) == hasseWittInvariant(f, p)) then (
		return false;
		)
	    else (return true);
	    )
	else (
	 -- now use the criteria for n=4
	-- need to check (d=1 and (-1,-1) != =hasseWittInvariant for f) for form to be anisotropic 
	    if (n==4) then (
		if ((equalUptoPadicSquare( d, 1, p)) and  (not (HilbertSymbol(-1,-1,p) == hasseWittInvariant(f, p)))) then (
		    return true;
		    )
		else (return false);
		)
	    else (
		print "Error: rank should have been 2, 3, 4, but isn't";
		return false;
		)
	    )
	)
    );

	    



