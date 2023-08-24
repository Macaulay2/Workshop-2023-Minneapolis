

-- Over k = R, the real symbol is 1 if either a or b is positive, and -1 if they are both negative.
-- See Serre, III Theorem 1

HilbertSymbolReal = method()
HilbertSymbolReal (QQ, QQ) := (ZZ) => (a, b) -> (
    if (a < 0 and b < 0) then (
	return -1;
	)
    else (
	return 1;
	);
    );

HilbertSymbolReal (QQ, ZZ) := (ZZ) => (a,b) -> (
    b1 := b/1;
    return HilbertSymbolReal(a, b1);
    );

HilbertSymbolReal (ZZ, QQ) := (ZZ) => (a,b) -> (
    a1 := a/1;
    return HilbertSymbolReal(a1, b);
    );

HilbertSymbolReal (ZZ, ZZ) := (ZZ) => (a,b) -> (
    a1:= a/1;
    b1:= b/1;
    return HilbertSymbolReal(a1, b1);
    );


-- Input: Any integers a and b and a prime p. The integers a and b are considered as elements of QQ_p.
-- Output: The Hilbert symbol (a,b)_p following Serre III Theorem 1

HilbertSymbol = method()
HilbertSymbol (ZZ, ZZ, ZZ) := (ZZ) => (a, b, p) -> (
    if (not isPrime(p)) then error "third argument of HilbertSymbol must be a prime";
    
    alpha := PadicValuation(a,p);
    beta := PadicValuation(b,p);
    u := sub(a/p^alpha,ZZ);
    v := sub(b/p^beta, ZZ);
    
    -- If epsilon(p) = 0
    if (p % 4 == 1) then(
	return ((squareSymbol(u,p))^beta * (squareSymbol(v,p))^alpha);
	);
    
    -- If epsilon(p) = 1
    if (p % 4 == 3) then(
	return ((-1)^(alpha*beta))*(((squareSymbol(u,p))^beta) * ((squareSymbol(v,p))^alpha));
	
	);
    
    -- Finally if p=2
    if p == 2 then(
	u = (u % 2);
	v = (v % 2);
	alpha = (alpha % 2);
	beta = (beta % 2);
	d := sub(((u-1)/2)*((v-1)/2) + alpha*((v^2-1)/8) + beta*((u^2-1)/8),ZZ);
	return ((-1)^d)
	);
    
    );


    
-- -- Given (a,b,p) integers, with a,b considered as p-adic numbers and p a prime, returns the Hilbert symbol (a,b)_p
-- HilbertSymbol = method()
-- HilbertSymbol (ZZ, ZZ, ZZ) := (ZZ) => (a, b, p) -> (
--     -- When p is odd, the Hilbert symbol (a,b)_p can be expressed by a formula using Legendre symbols
--     -- or equivalently, be defined, by the following condition.  
--     if (odd p) then (
-- 	if (squareSymbol(a, p)==1 or squareSymbol(b, p)==1 or 
-- 	    squareSymbol(-1*a*b, p) ==1 or
-- 	    (squareSymbol(a,p)*squareSymbol(b,p)==1)) then (
-- 		hilb:=1;
-- 		)
-- 	else (
-- 	    hilb=-1;
-- 	    );
-- 	)
--      else (
-- 	 e1:=PadicValuation(a, p);
-- 	 u:=sub(a/p^e1,ZZ);
-- 	 e2:=PadicValuation(b, p);
-- 	 b1:=sub(b/p^e2, ZZ);
-- 	 c1:= (u-1)/2;
-- 	 c2:= (b1-1)/2;
-- 	 d1:= (u^2-1)/8;
-- 	 d2:= (b1^2-1)/8;
	 
	 
-- 	 d:= c1*c2+e1*d2+e2*d1;
	  
-- 	 -- when p=2, the Hilbert symbol (a,b)_2 equals (-1)^d
-- 	 if (even sub(d,ZZ)) then (
-- 	     return 1;
-- 	     )
-- 	 else (
-- 	     return -1;
-- 	     );
-- 	 	    );
--     return hilb;	
-- );

HilbertSymbol (QQ, QQ, ZZ) := (ZZ) => (a, b, p) -> (
 
-- if a, b are rational numbers with denominators,  one can multiply by square of denominator to 
-- get a', b' integers.  Then evaluate hilbertSymbol (a', b', p);
    
    if not liftable(a, ZZ) then (
	a1 := numerator a;
	a2 := denominator a;
	a = a1*a2;
	);
	
    if not liftable(b, ZZ) then (
	b1 := numerator b;
	b2 := denominator b;
	b = b1*b2;
	);
    
    a = sub(a,ZZ);
    b = sub(b,ZZ);
    return HilbertSymbol(a,b,p);
    );

HilbertSymbol (ZZ, QQ, ZZ) := (ZZ) => (a, b, p) -> (
    a1:=a/1;
    return HilbertSymbol(a1,b, p);    
   );

HilbertSymbol (QQ, ZZ, ZZ) := (ZZ) => (a, b, p) -> (
   b1:=b/1;
   return HilbertSymbol(a, b1, p);    
   );


