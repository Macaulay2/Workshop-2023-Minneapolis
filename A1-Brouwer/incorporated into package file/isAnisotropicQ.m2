--Nikita, Tom, Zhaobo

path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");

needs "GW-type.m2"
needs "getInvariants.m2"
needs "diagonalize.m2"
needs "squarefreepart.m2"
needs "discForm.m2"
needs "isAnisotropicDiagQp.m2"



--isAnisotropicQ takes n GWClass over QQ and returns a Boolean based on whether or not 
--the class is anisotropic
--unlike simplifyForm it can say for certain if the form is anisotropic or not since it
--does not use rationalPoints

-- Input:  A GrothendieckWittClass for a quadratic form
-- Output: True if form is Anisotropic;  False, if form is Isotropic


isAnisotropicQ = method()
isAnisotropicQ (GrothendieckWittClass) := Boolean => (alpha) -> (
    A:= alpha.matrix;
    n:= numRows(A);
    kk:= ring A;

    if (not (kk===QQ)) then (error "GrothendieckWittClass is not over QQ");

    --check if form is degenerate
    signature := (getInvariants(alpha))_2;
    if (signature_1 > 0) then (return false);
    
     -- if rank =1, then a non-degenerate form is anisotropic
    if (n==1) then (return true);
    
    --if rank>=5, we can use signature do decide this
    --the non-degenerate form will be anisotropic iff all diagonal entries have same sign
    if (n>= 5) then (
        return ((signature_0 == n) or (signature_2 == n));
    );
   

    --if 2<= rank <=4, we need to take p-adic completions
    -- First, we diagonalize the matrix
    diagA := diagonalize(A);  
    -- Then obtain the diagonal entries.  These will be rational numbers;
    diagEntriesA := apply(n, i-> A_(i,i));
    -- Make then integers by multiplying by integer squares (so that the forms are equivalent);
    -- The sub command forces the list to be integers
    diagIntEntriesA:= apply(n, i-> squarefreePart(sub(numerator(diagEntriesA_i) * denominator(diagEntriesA_i),ZZ)));
    
    -- disc = discriminant of form, product of diagonal elements.  
    disc = discForm(diagIntEntriesA);
    
    
   
   
    -- Using Q_p criteria from Thm 6, Section 2.2 of Serre's Course in Arithmetic
   
    -- For n=2, need -disc to be a square for all p to be isotropic, so in particular, need -disc=1 for isotropic
    if (n==2) then (  
	 -- Make disc a squarefree integer
	 d2 := squarefreePart(disc); 
	 return (not (d2==-1) ) 
	 );
     
  
    -- if p>2, then hilbert symbol (a,b)=1 if, a, b not divisible by p.  So hasseWittInvariant is also 1.  
    -- Then for n=3,4, form is automatically isotropic if p doesn't divide disc. 
    -- So only need to check if f is anisotropic over Q_p for p=2 and primes p dividing disc.
    
    -- first check p=2 case
    if (isAnisotropicDiagQp(diagIntEntriesA, 2)) then (return true);
    
    -- create list of primes dividing disc
    -- first take absolute value of disc
    d1:= disc;
    if (d1<0) then (d1=-d1);
    
    -- H is HashTable of factors of disc
    
    H:= hashTable( factor d1);
    -- the keys k are the prime factors
    k:= keys H;
    i:=0;
   
    while (i< #k) do (
	    	p = k_i;
		if (isAnisotropicDiagQp(diagIntEntriesA, p)) then (return true);
		i=i+1;
		);

-- if the function hasn't returned false yet, then isotropic over all primes p.  hence form is isotropic over Q	    
  
    return false;
);
	
	
