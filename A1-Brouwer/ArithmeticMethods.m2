
------------------------
-- Arithmetic operations
------------------------

-- Input: A rational number
-- Output: Smallest magnitude integer in its square class

squarefreePart = method()
squarefreePart (QQ) := (ZZ) => (n) -> (
    if n == 0 then (
        return 0
        );
    if n > 0 then (
        tableOfPrimeFactorsQQ:=hashTable(factor(numerator(n)*denominator(n)));
        return product(apply(keys(tableOfPrimeFactorsQQ),p->p^(tableOfPrimeFactorsQQ#p%2)))
        );
    if n < 0 then (
        tableOfPrimeFactorsQQNeg:=hashTable(factor(numerator(-n)*denominator(-n)));
        return -product(apply(keys(tableOfPrimeFactorsQQNeg),p->p^(tableOfPrimeFactorsQQNeg#p%2)))
        );
    )

squarefreePart (ZZ) := (ZZ) => (n) -> (
    if n == 0 then (
        return 0
        );
    if n > 0 then (
        tableOfPrimeFactors:=hashTable(factor(n));
        return product(apply(keys(tableOfPrimeFactors),p->p^(tableOfPrimeFactors#p%2)))
        );
    if n < 0 then (
        tableOfPrimeFactorsNeg:=hashTable(factor(-n));
        return -product(apply(keys(tableOfPrimeFactorsNeg),p->p^(tableOfPrimeFactorsNeg#p%2)))
        );
    )

-- Returns a list of prime factors of an integer
primeFactors = method()
primeFactors (ZZ) := List => (n) -> (
    if abs(n)==1 then(
	return {};
	);
    
    return keys(hashTable(factor(abs(n))))
    );

primeFactors (QQ) := List => (n) -> (
    if (not liftable(n,ZZ) == true) then(
	error "tried to take prime factors of a rational";
	);
    
    return primeFactors(sub(n,ZZ));
    
    )


-- Input: Element of a finite field
-- Output: True if an element of a finite field is a square, false otherwise
legendreBoolean = method()
legendreBoolean (RingElement) := (Boolean) => a -> (
    if not instance(ring(a),GaloisField) then error "Error: this works only for Galois fields";
    q := (ring(a)).order;
    -- Detects if a is a square in F_q
    a^((q-1)//2) == 1 
    )

------------------------------
-- Commutative algebra methods
------------------------------

-- Input: A list L of functions f1,...,fn over the same ring R and p is a prime ideal of an isolated zero.
-- Output: A list of basis elements of the local k-algebra Q_p(f) = R[x]_p/(f).

localAlgebraBasis = method()
localAlgebraBasis (List, Ideal) := (List) => (L,p) -> (
    
    -- Determine whether or not an ideal is prime
    if isPrime(p) == false then (
        error "Error: ideal is not prime"
        );
    
    -- Ambient ring
    R := ring L#0;
    I := ideal(L);
    
    -- Check whether or not an ideal is zero-dimensional
    if dim I > 0  then (
        error "Error: morphism does not have isolated zeroes"
        );
    if (not isSubset(I,p)) then (
        error "Error: prime is not a zero of function"
        );
    J := I:saturate(I,p);
    A := R/J;
    B := basis(A);
    return flatten(entries(B))
    )

-- Input: A zero-dimensional ideal (f_1,..f_n) < k[x_1..x_n].
-- Output: The rank of the global algebra  k[x_1..x_n].

rankGlobalAlgebra = method()
rankGlobalAlgebra (List) := (ZZ) => (Endo) -> (
    
    -- Get the underlying field    
    kk := coefficientRing(ring(Endo#0));    
    if isField(kk) == false then(
    	kk = toField(kk);
    	);
    
    -- Let S = k[x_1..x_n] be the ambient polynomial ring
    S := ring(Endo#0);
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then (
	error "Error: ideal is not zero-dimensional";
	);
    
    -- Get the rank of S/ideal(Endo) as a kk-vector space
    return numColumns(basis(S/ideal(Endo)));   
    )



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
	x := getSymbol "x";
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

-- Check if something is a p-adic square
isPadicSquare = method()
isPadicSquare (ZZ, ZZ):= (Boolean) => (a, p) -> (
    return equalUptoPadicSquare(a,1,p)
    );
