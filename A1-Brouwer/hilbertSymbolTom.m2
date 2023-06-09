-- Calculate the Hilbert symbol in Q_p

--  Very preliminary version of Invariant Calculations (hilbert symbol, etc.) for quadratic functions over Q
-- 6/9/23
-- Tom Hagedorn


-- Input: Integers a, b, and prime p 
-- Output: The Hilbert symbol, (a,b)_p

load "simplifyForm.m2"

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

-- Given a number a and a prime p
-- Output: 1, if a is a unit square,  -1, if a=p^even (unit square), 0 otherwise
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
	 c1:= eps2Hilbert(a1);
	 c2:= eps2Hilbert(b1);
	 d:= c1*c2+e1*c2+e2*c1;
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

    
 -- Two Q forms over Q_p are isomorphic if they have same rank, same discriminant, and same eps
 -- invariant   
    
    
    
      
epsilonHilbert = method()

epsilonHilbert(list) := ZZ => f -> (
       a:=1;
       len:=#f
       for i from 0 to len -1 do (
	   if (not ring(f_i) == ZZ) then print"Error:  Hilbert symbol evaluated at a non-integer";
       for i from 0 to len-2 do (
       	   for j from i+1 to len-1 do (
	       a= a* hilbertSymbol(f_i, f_j);
	       );
	   );
       return a;          
    );
    
