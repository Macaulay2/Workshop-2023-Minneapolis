
-- discForm computes the product of the entries of a list;
-- used to compute the discriminant of a diagonal form, where form is given by list of diagonal entries
discForm = method ()
discForm (List):=(ZZ) => (f) -> (
    -- Input: f = list of diagonal elements of quadratic form
    -- Output: Product of the entries
    disc:=1;
    --for loop counts the number of negative diagonal entries
    for i from 0 to (#f -1) do(
	 disc = disc* f_i;
	 );
    return disc;
    );

discForm (Sequence):=(ZZ) => (f) -> (
    -- Input: f = list of diagonal elements of quadratic form
    -- Output: Product of the entries
    disc:=1;
    --for loop counts the number of negative diagonal entries
    for i from 0 to (#f -1) do(
	 disc = disc* f_i;
	 );
    return disc;
    );

-- Calculate the Hilbert symbol in Q_p

--  Very preliminary version of Invariant Calculations (hilbert symbol, etc.) for quadratic functions over Q
-- 6/9/23
-- Tom Hagedorn
-- Currently, all calculations are done for a  quadratic form, diagonalized, and we work with an isomorphic rational form, scaled 
-- so that all diagonal elements are integers.

-- Goals  (6/9/23)
-- Done:
-- 0. Expand code to have function that outputs all the invariants for a rational quadratic form over Qp 
-- 1.  Use invatiants of two forms to tell whether two rational forms are isomorphic over Q
-- 2.  Use invariant to determine if a rational q. form is anisotropic

-- Todo:
-- 1.  Expand code to have function that outputs all the invariants for a rational quadratic form.  
-- 3.  Expand code so that it works for any quadratic form over a number field (should be very similar to the code below).


-- Input: Integers a, b, and prime p 
-- Output: The Hilbert symbol, (a,b)_p

exponentPrimeFact = method()
exponentPrimeFact (ZZ, ZZ) := (ZZ) => (n, p) -> (
    if (n<0) then (n=-n);
    if (n==0) then error "Error: Trying to find prime factorization of 0";
    H:=hashTable (factor n);
    a:=0;
    if H#?p then (
    	a=H#p;)
    else (
	a=0;
	);
    return a;
    );




-- Input: An integer a and a prime p
-- Output: 1, if a is a unit square,  -1, if a=p^(even power)x  (non-square unit), 0 otherwise
-- Note:  The terminology "Square Symbol" comes from John Voight's Quaternion Algebra book
squareSymbol = method()
squareSymbol(ZZ, ZZ) := (ZZ) => (a, p) -> (
    x := getSymbol "x";
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
    
-- Given (a,b,p) integers, with a,b considered as p-adic numbers and p a prime, returns the Hilbert symbol (a,b)_p
HilbertSymbol = method()
HilbertSymbol (ZZ, ZZ, ZZ) := (ZZ) => (a, b, p) -> (
    -- When p is odd, the Hilbert symbol (a,b)_p can be expressed by a formula using Legendre symbols
    -- or equivalently, be defined, by the following condition.  
    if (odd p) then (
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
	 c1:= (a1-1)/2;
	 c2:= (b1-1)/2;
	 d1:= (a1^2-1)/8;
	 d2:= (b1^2-1)/8;
	 
	 
	 d:= c1*c2+e1*d2+e2*d1;
	  
	 -- when p=2, the Hilbert symbol (a,b)_2 equals (-1)^d
	 if (even sub(d,ZZ)) then (
	     return 1;
	     )
	 else (
	     return -1;
	     );
	 
	 --need to do even case
	    );
    return hilb;	
);

HilbertSymbol (QQ, QQ, ZZ) := (ZZ) => (a, b, p) -> (
    if not liftable(a, ZZ) then error "first argument must be an integer";
    if not liftable(b, ZZ) then error "second argument must be an integer";
    a = sub(a,ZZ);
    b = sub(b,ZZ);
    return HilbertSymbol(a,b,p)	
)


HilbertSymbolReal = method()
HilbertSymbolReal (ZZ, ZZ):=(ZZ) => (a,b)->(
    if (a<0 and b<0) then (
	return -1;
	)
    else (
	return 1;
	)
    );
