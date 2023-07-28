

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

-- Calculate the Hilbert symbol in Q_p

--  Very preliminary version of Invariant Calculations (hilbert symbol, etc.) for quadratic functions over Q
-- 6/9/23
-- Tom Hagedorn
-- Currently, all calculations are done for a  quadratic form, diagonalized, and we work with an isomorphic rational form, scaled 
-- so that all diagonal elements are integers.

-- Goals  (6/9/23)
-- Done:
-- 0. Expand code to have function that outputs all the invariants for a rational quadratic form over Qp 
--     Done with function invariantFormQp
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


signatureRealQForm = method ()

signatureRealQForm (List):=(ZZ, ZZ, ZZ) => (f) -> (
    -- Input: f = list of diagonal elements of quadratic form
    -- Output: (r, s, t): Signature of form over the reals. 
    --      r= # of positive entries, s= # of negative entries, t= number of zeros
    posEntries :=0;
    negEntries:= 0;
    zeroEntries:=0;
    --for loop counts the number of negative diagonal entries
    n:=#f;
    for i from 0 to (n-1) do(
            if f_(i)>0 then(
                posEntries=posEntries+1;
            	)
	    else (
		if f_(i)<0 then(
                    negEntries=negEntries+1;
            	)
	    else (
		zeroEntries = zeroEntries +1;
		)
	    )
	);
    return (posEntries, negEntries, zeroEntries);
    
	     
    );





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
	    if (d==0) then print"Error: Form is degenerate";
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


isIsomorphicFormQ = method ()

isIsomorphicFormQ (Matrix, Matrix):= (Boolean) => (f, g) -> (
    -- Input: (f,g):  f, g are two square bilinar forms over Q, represented by matrices
    -- Output: true if f, g are isomorphic over Q
    
    -- Check same size
    if (numRows f != numRows g) then (return false;);
    
    -- First, we diagonalize both matrices
    df:= diagonalize f;
    dg:= diagonalize g;
    
   
    -- Then make all entries integers, by multiplying by denominator squared, and clearing squares
    -- create list of diagonal entries in integers
    n:=numRows f;
    f1:=apply(n, i-> df_(i,i));
    g1:=apply(n, i-> dg_(i,i));
    
  
    -- convert rational diagonals to square-free integers by multiplying by squares
    
    f2:= apply(n, i-> squarefreePart(sub(numerator(f1_i) * denominator(f1_i),ZZ)));
    g2:= apply(n, i-> squarefreePart(sub(numerator(g1_i) * denominator(g1_i),ZZ)));
    
    print f2;
    print g2; 
    -- Now compare forms
    return isIsomorphicDiagFormQ(f2, g2);
    
    );
    
    
    
isIsotropicDiagFormQp = method()

isIsotropicDiagFormQp (List, ZZ):=(Boolean) => (f, p) -> (
    -- Input: (f, p): Diagonal Form with Integer coeffs, p a prime
    -- Output: true if form represents 0 over Qp
    
    n:=#f;
    d:=discForm(f);
    e:=hasseWittInvariant(f, p);
    if (n==2) then (
    -- need to compare d ==-1 in Qp^*/Qp^2	
	if (equalUptoPadicSquare(sub(-1, ZZ), d, p)) then (
	    return true;
	    )
	else (return false);
	)
    else (
	if (n==3) then (
	    if (hilbertSymbol(-1, -d, p) == hasseWittInvariant(f, p)) then (
		return true;
		)
	    else (return false);
	    )
	else (
	    if (n==4) then (
		if ((not equalUptoPadicSquare( d, 1, p)) or (equalUptoPadicSquare(d, 1, p) and (hilbertSymbol(-1,-1,p) == hasseWittInvariant(f, p)))) then (
		    return true;
		    )
		else (return false);
		)
	    else (
		return true;
		)
	    )
	)
    );


isIsotropicDiagFormQ = method()


isIsotropicDiagFormQ (List):=(Boolean) => (f) ->(
    -- need to check if anisotropic over RR and all Qp
    n:= #f;
    a:=signatureRealQForm(f);
    -- if real form is definite then form f is not isotropic
    if (a_0 * a_1 <=0 ) then (
	return false;
	) 
    else (
	if (n>4) then (return true; ) 
	else (
	    if (n==1 ) then (return false; );
	    -- for n=2,3,4, the form is isotropic if p!=2 and p doesn't divide any of the 
	    -- diagonal terms as the conditions in Q_p are automatically satisfied.
	    d:= discForm(f);
	    d1:= d;
	    if (d1<0) then (d1=-d1);
	    H:= hashTable( factor d1);
	    k:= keys H;
	    i:=0;
	    l:= #k;
	    if (not isIsotropicDiagFormQp(f, 2)) then (return false);
	    while (i< l) do (
	    	p = k_i;
		if ( not isIsotropicDiagFormQp(f, p)) then (return false);
		i=i+1;
		);
	    return true;
	    );
	);
    );	
			
			
	
	
	
invariantFormQ =method()

-- Input: (Q): Quadratic form Q given by list of diagonal elements.  
--  	Assume list constists of integers
-- Output:  (Rank, Disc, Signature, Hasse Invariant for all primes p when not 1)

invariantFormQ (List):= (ZZ, ZZ, List, List) => (f) -> (
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
    b:=discForm(f);
    c:=signatureRealQForm(f);
    d:=b;
    if (b<0) then (d=-b);
    -- The keys of H contain all primes dividing coefficients
    H:= hashTable( factor d);
    k:= keys H;
    l:={};
    if ((not H#?2) and hasseWittInvariant(f, 2) == -1) then l=append(l, 2);
    for i from 0 to #k-1 do (
	if  (hasseWittInvariant(f, k_i) == -1) then (
	    l=append(l,k_i);
	    );
	);
    l=sort(l);
    return (a, b, c, l);
    
    );
  
