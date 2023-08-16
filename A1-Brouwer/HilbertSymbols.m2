
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


-- Computes the Hilbert symbol (a,b)_p following Serre III Theorem 1
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
    if p==2 then(
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
    if not liftable(a, ZZ) then error "first argument must be an integer";
    if not liftable(b, ZZ) then error "second argument must be an integer";
    a = sub(a,ZZ);
    b = sub(b,ZZ);
    return HilbertSymbol(a,b,p)	
)

-- Over k=R, the real symbol is 1 if either a or b is positive, and -1 if they are both negative.
-- See Serre, III Theorem 1
HilbertSymbolReal = method()
HilbertSymbolReal (ZZ, ZZ):=(ZZ) => (a,b)->(
    if (a<0 and b<0) then (
	return -1;
	)
    else (
	return 1;
	)
    );
