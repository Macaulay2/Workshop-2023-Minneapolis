---------------------------------------
-- Invariants
---------------------------------------

-- Input: A Grothendieck-Witt class beta defined over QQ or RR
-- Output: The number of positive entries in its diagonalization

numPosEntries = method()
numPosEntries (GrothendieckWittClass) := ZZ => beta ->(
    B := beta.matrix;
    n := numRows(B);
    kk := ring B;
    if not (kk === RR or instance(kk,RealField) or kk === QQ) then(
        error "Field is not QQ or RR";
        );
    diagB := congruenceDiagonalize(B);
    posEntries := 0;
    for i from 0 to (numRows(B) - 1) do (
        if diagB_(i,i) > 0 then(
            posEntries = posEntries+1;
            );
	);

    return posEntries
    );

-- Input: A Grothendieck-Witt class beta defined over QQ or RR
-- Output: The number of negative entries in its diagonalization

numNegEntries = method()
numNegEntries (GrothendieckWittClass) := ZZ => beta ->(
    B := beta.matrix;
    n := numRows(B);
    kk := ring B;
    if not (kk === RR or instance(kk,RealField) or kk === QQ) then(
        error "Field is not QQ or RR";
        );
    diagB := congruenceDiagonalize(B);
    negEntries := 0;
    for i from 0 to (numRows(B)-1) do (
        if diagB_(i,i) < 0 then(
            negEntries = negEntries+1;
            );
	);

    return negEntries
    );

-- Input: A Grothendieck-Witt class beta defined over QQ or RR
-- Output: The signature of beta

signature = method()
signature (GrothendieckWittClass) := ZZ => (beta) ->(
    sig := numPosEntries(beta) - numNegEntries(beta);
    return sig
    );

---------------------------
-- Comparing forms over QQ
---------------------------

-- Input: A Grothendieck-Witt class beta defined over QQ
-- Output: A squarefree integral representative of its discriminant

integralDiscriminant = method()
integralDiscriminant (GrothendieckWittClass) := (ZZ) => (beta) -> (
    B:= beta.matrix;
    rankForm:= numRows(B);
    kk:= ring B;
    
    if (not (kk === QQ)) then (error "GrothendieckWittClass is not over QQ");
    
    -- Take an integral diagonal representative for beta
    gamma := integralDiagonalRep(beta);
    G := gamma.matrix;
    
    discrimForm:= 1;
    for i from 0 to (rankForm-1) do(
	discrimForm = discrimForm * (G_(i,i));
	);
    
    return sub(squarefreePart(discrimForm),ZZ);
    );

-- Input: A Grothendieck-Witt class beta defined over QQ
-- Output: The smallest list of primes that divide its discriminant

relevantPrimes = method()
relevantPrimes (GrothendieckWittClass) := List => (beta) -> (
    B := beta.matrix;
    rankForm := numRows(B);
    kk:= ring B;
    
    -- Take a diagonal integral representative of the form
    gamma := integralDiagonalRep(beta);
    D := diagonalEntries(gamma);
    
    L := {};
    
    -- Append all the prime factors of each of the entries appearing on a diagonal
    for x in D do(
	L = unique(L | primeFactors(sub(x,ZZ)));
	);
    
    return L
    );

-- Two Q forms over Q_p are isomorphic if they have same rank, same discriminant, and same Hasse-Witt invariant   

HasseWittInvariant = method()

-- epsilonHilbert computes the epsilon function for a diagonal quadratic form over Q_p
-- Function requires the list of the diagonal elements of the quadratic form, to be integers

-- Input:  A list of the diagonal elements (f_i) for the quadratic form, assumed to be integers, and a prime p
-- Output: The HasseWittInvariant function for the quadratic form (f_i) for Q_p

HasseWittInvariant (List, ZZ) := ZZ => (L,p) -> (
       a := 1;
       len := #L;
       
       -- Replace every entry of L with its squarefree part so we can be sure we're evaluating at integers
       f := {};
       for x in L do(
	   f = append(f,squarefreePart(x));
	   );
       
       for i from 0 to len - 1 do (
	   if not liftable(f_i,ZZ) then (error "Error:  Hilbert symbol evaluated at a non-integer");
	   );
       for i from 0 to len - 2 do (
       	   for j from i + 1 to len - 1 do (
	       a = a * HilbertSymbol(f_i, f_j, p);
	       );
	   );
       
       return a;
    );

HasseWittInvariant(GrothendieckWittClass, ZZ) := ZZ => (beta,p) -> (
    kk := baseField beta;
    if not (kk === QQ) then error "method is only implemented over the rationals";
    if not isPrime(p) then error "second argument must be a prime number";
    return HasseWittInvariant(diagonalEntries(beta),p)
    
    )
