-- Calculate the Hilbert symbol in Q_p

--  Very preliminary version of Invariant Calculations (hilbert symbol, etc.) for quadratic functions over Q
-- 6/9/23
-- Tom Hagedorn
-- Currently, all calculations are done for a quadratic form, diagonalized, and we work with an isomorphic rational form, scaled 
-- so that all diagonal elements are integers.

-- Goals  (6/9/23)
-- 0. Expand code to have function that outputs all the invariants for a rational quadratic form over Qp 
--     Done with function invariantFormQp
-- 1.  Expand code to have function that outputs all the invariants for a rational quadratic form.  
-- 2.  Use invariant to determine if a rational q. form is anisotropic
-- 3.  Use invatiants of two forms to tell whether two rational forms are isomorphic over Q
-- 4.  Expand code so that it works for any quadratic form over a number field (should be very similar to the code below).


-- Input: Integers a, b, and prime p 
-- Output: The Hilbert symbol, (a,b)_p

load "simplifyForm.m2"
load "squarefreepart.m2"

exponentPrimeFact = method()

exponentPrimeFact (ZZ, ZZ) := (ZZ) => (n, p) -> (
    if (n<0) then (n=-n);
    if (n==0) then print"Error: Trying to find prime factorization of 0";
    H:=hashTable (factor n);
    if H#?p then (
    	a:=H#p;)
    else (
	a=0;
	);
    return a;
    );


squareSymbol = method()

-- Input: An integer a and a prime p
-- Output: 1, if a is a unit square,  -1, if a=p^(even power)x  (non-square unit), 0 otherwise
-- Note:  The terminology "Square Symbol" comes from John Voight's Quaternion Algebra book

squareSymbol(ZZ, ZZ) := (ZZ) => (a, p) -> (
    R:=GF(p, Variable => x);
    e1:=exponentPrimeFact(a,p);
    if (even e1) then (
    	a1:= sub(a/(p^e1), ZZ);
	a2:=sub(a1, R);
	if legendreBoolean(a2) then (
	    ans:=1;
	    ) 
	else (
	    ans=-1;
	    );
	)
    else (
	ans=0;
	);
    return ans;
    );
    
    
eps2Hilbert = method()

eps2Hilbert(ZZ):= ZZ => (a) -> (
	    if (odd a) then (
		return (a-1)/2;
	        )
	    else (
		print "Error: epsilonHilbert applied to even integer";
	    );    
    );

omegaHilbert = method()

omegaHilbert(ZZ):= ZZ  => (a) -> (
    if (odd a) then (
		return (a^2-1)/8;
	        )
	    else (
		print "Error: omegaHilbert applied to even integer";
	    );    
    
    );

hilbertSymbol = method()

hilbertSymbol (ZZ, ZZ, ZZ) := (ZZ) => (a, b, p) -> (
-- Note: The Hilbert symbol (a,b)_p is defined for a, b p-adic numbers
-- We will be assuming that the a, b are integers.
-- We need to distinguish 
    
    if (odd p) then (
	-- when p is odd, the Hilbert symbol (a,b)_p can be expressed by a formula using Legendre symbols
	-- or equivalently, be defined, by the follwoing condition.  
	if (squareSymbol(a, p)==1 or squareSymbol(b, p)==1 or 
	    squareSymbol(-1*a*b, p) ==1 or
	    (squareSymbol(a,p)*squareSymbol(b,p)==1)) then (
		hilb:=1;
		)
	else (
	    hilb=-1;
	    );
	)
     else (
	 e1:=exponentPrimeFact(a, p);
	 a1:=sub(a/p^e1,ZZ);
	 e2:=exponentPrimeFact(b, p);
	 b1:=sub(b/p^e2, ZZ);
	 c1:= eps2Hilbert(a1);
	 c2:= eps2Hilbert(b1);
	 d1:= omegaHilbert(a1);
	 d2:= omegaHilbert(b1);
	 
	 d:= c1*c2+e1*d2+e2*d1;
	 -- when p=2, the Hilbert symbol (a,b)_2 equals (-1)^d
	 if (even d) then (
	     return 1;
	     )
	 else (
	     return -1;
	     );
	 
	 --need to do even case
	    );
    return hilb;	
);

hilbertSymbolReal = method()

hilbertSymbolReal (ZZ, ZZ):=(ZZ) => (a,b)->(
    if (a<0 and b<0) then (
	return -1;
	)
    else (
	return 1;
	)
    );

    
 -- Two Q forms over Q_p are isomorphic if they have same rank, same discriminant, and same Hasse-Witt invariant   
    
    
    
      
hasseWittInvariant = method()

-- epsilonHilbert computes the epsilon function for a diagonal quadratic form over Q_p
-- Function requires the list of the diagonal elements of the quadratic form, to be integers
-- Input:  A list of the diagonal elements (f_i) for the quadratic form, assumed to be integers, and a prime p
-- Output: The hasseWittInvariant function for the quadratic form (f_i) for Q_p

hasseWittInvariant (List, ZZ) := ZZ => (f,p) -> (
       a:=1;
       len:=#f;
       for i from 0 to len-1 do (
	   if (not ring(f_i) === ZZ) then (print"Error:  Hilbert symbol evaluated at a non-integer");
	   );
       for i from 0 to len-2 do (
       	   for j from i+1 to len-1 do (
	       a= a * hilbertSymbol(f_i, f_j, p);
	       );
	   );
       
       return a;          
    );


invariantFormQp =method()

-- Input: (Q, p): Quadratic form Q given by list of diagonal elements.  For now, assume list to to be integers
--                p is a prime number
-- Output:  (Rank, Disc, Hasse Invariant)

invariantFormQp (List, ZZ):= (ZZ, ZZ, ZZ) => (f, p) -> (
    -- currently will export the discriminant as a square free integer
    -- Note:  Still need a way to treat two integers as defining the same discriminant, if they differ by a
    --  square in Z_p.  They need to have the same parity of power of prime p, and the quotient of their 
    -- (prime-to-p) parts must define a unit square in Z_p^*.

    len:=#f;
    for i from 0 to (len-1) do (
	if (f_i==0) then (print"Error: Form is degenerate");
	if (not ring(f_i)===ZZ ) then (print "Error: Diagonal elements of form should be integers");
	);
    a:=len;
    b:=1;
    for i from 0 to len-1 do (b=b*f_i);
    b=squarefreePart(sub(b, QQ));
    c:=hasseWittInvariant(f, p);
    return(a, b, c);
    );

-- Need routine 

equalUptoPadicSquare = method()

equalUptoPadicSquare (ZZ, ZZ, ZZ):= (Boolean) => (a, b, p) -> (
-- Given a, b integers, determines if a, b differ by a square in Q_p
-- One has to handle the cases when p is odd, and p=2 differently
if (odd p) then (
    -- p is odd and we need to check that the powers of p have the same parity, and the units
    -- differ by a square in GF(p)
    a1:=squarefreePart(a);
    b1:=squarefreePart(b);
    if (exponentPrimeFact(a1, p ) != exponentPrimeFact(b1, p)) then (
	return false;
        )
    else (
    	-- c1 will be an integer prime to p
	c1:= squarefreePart(a1*b1);
	return (legendreBoolean( sub(c1, GF(p, Variable => x)))); 
	);
    )
else (
    -- Case when p=2.  Then we have to check that the powers of p have the same parity, and 
    -- that the units agree mod 8.
    a1=squarefreePart(a);
    b1=squarefreePart(b);
    if (exponentPrimeFact(a1, p ) != exponentPrimeFact(b1, p)) then (
	return false;
        )
    else (
    	-- c1 will be an integer prime to p
	c1= squarefreePart(a1*b1);
	c1 = c1 % 8;
	-- if c1 =1, then the two odd units are congruent mod 8, and are squares in Q2
	return (c1==1); 
	);
    );
  );  
 
isIsomorphicFormQp = method ()

isIsomorphicFormQp (List, List, ZZ):= (Boolean) => (f, g, p) -> (
    -- Input: (f,g ,p):  f, g are lists of diagonal elements (in integers) of two forms over Qp
    --         p is a prime number
    -- Output: true if f, g are isomorphic over Qp
    
    -- Should check that the elements in the list are all integers. 
    a:= invariantFormQp(f,p);
    b:= invariantFormQp(g,p);
    if (a_0 == b_0 and a_2 == b_2 and  equalUptoPadicSquare(a_1, b_1, p)) then (
	-- all three of the invariant are the same, so the forms are isom. over Qp
	return true;
	)
    else (
	-- one of the invariants differ, so the forms are different.
	return false;
    );
);
