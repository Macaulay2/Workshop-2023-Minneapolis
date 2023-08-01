--A1BrouwerDegrees.m2
newPackage(
    "A1BrouwerDegrees",
    Version=>"1.0",
    Date=>"June 5, 2023",
    Authors=>{
        {Name=>"Nikita Borisov",
	 Email=>"nborisov@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/nikita-borisov"},
        {Name=>"Thomas Brazelton",
	 Email=>"tbraz@math.upenn.edu",
	 HomePage=>"https://www2.math.upenn.edu/~tbraz/"},
        {Name=>"Frenly Espino",
	 Email=>"frenly@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/frenly-espino"},
         {Name=>"Tom Hagedorn",
	 Email=>"hagedorn@tcnj.edu",
	 HomePage=>"https://hagedorn.pages.tcnj.edu/"},
        {Name=>"Zhaobo Han",
	 Email=>"zbtomhan@sas.upenn.edu",
	 HomePage=>"https://www.linkedin.com/in/zhaobo-han-77b1301a2/"},
     	{Name=>"Jordy Lopez Garcia",
	 Email=>"jordy.lopez@tamu.edu",
	 HomePage=>"https://jordylopez27.github.io/"},
        {Name=>"Joel Louwsma",
	 Email=>"jlouwsma@niagara.edu",
	 HomePage=>"https://www.joellouwsma.com/"},
        {Name=>"Andrew Tawfeek",
	 Email=>"atawfeek@uw.edu",
	 HomePage=>"https://www.atawfeek.com/"},
        {Name=>"Wern Juin Gabriel Ong",
	 Email=>"gong@bowdoin.edu",
	 HomePage=>"https://wgabrielong.github.io/"}
	},
    Headline=>"TODO",
    PackageImports=>{
	"RationalPoints2",
	"RealRoots"
	},
    PackageExports=>{},
    DebuggingMode=>true
    )


export{
    -- types
    "GrothendieckWittClass",
    -- arithmetic methods
    "squarefreePart",
    "legendreBoolean",
    -- matrix methods
    -- methods
    "baseField",
    "gwClass",
    "gwAdd",
    "gwMultiply",
    "diagonalize",
    "diagonalForm",
    "simplifyForm",
    "simplifyFormString",
    "globalA1Degree",
    "localA1Degree",
    "localAlgebraBasis",
    "signature",
    "isIsomorphic2",
    "hilbertSymbol",
    "isIsotropic",
    "isAnisotropic"
    }

--------------
-- Arithmetic operations
--------------

-- Takes in a rational number or integer and outputs the smallest magnitude integer in its square class
squarefreePart = method()
squarefreePart (QQ) := (ZZ) => (n) -> (
    if n==0 then (
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
    if n==0 then (
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

-- Boolean function, inputting an element of a finite field. Returns true if an element of a finite field is a square and false otherwise
legendreBoolean = method()
legendreBoolean (RingElement) := (Boolean) => a -> (
    if not instance(ring(a),GaloisField) then error "Error: this works only for Galois fields";
    q := (ring(a)).order;
    -- The following Boolean detects if a is a square in F_q
    a^((q-1)//2) == 1 
)

-------------------------------------------
-- Algebra methods
-------------------------------------------

-- Input: L is list of functions {f1,...,fn} over same ring R and p is prime ideal of an isolated zero

-- Output: list of basis elements of local k-algebra Q_p(f) where f = (f1,...,fn):A^n --> A^n

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

-- This is a method for obtaining the rank of a global algebra k[x_1..x_n] over a field, given the input of a zero-dimensional ideal (f_1,..f_n) < k[x_1..x_n]
rankGlobalAlgebra = method()
rankGlobalAlgebra (List) := (ZZ) => (Endo) -> (
    -- Get the underlying field    
    kk := coefficientRing(ring(Endo#0));    
    if isField(kk) == false then(
    	kk = toField(kk);
    	);
    
    -- Let S = k[x_1..x_n] be the ambient polynomial ring
    S:=ring(Endo#0);
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then (
	print "Error: ideal is not zero-dimensional";
	return Endo;
	);
    
    -- Get the rank of S/ideal(Endo) as a kk-vector space
    return numColumns(basis(S/ideal(Endo)));
    
    );



---------------
-- Matrix manipulations
---------------

-- Block sum "++" of a zero matrix with something else outputs the wrong thing
safeBlockSum = method()
safeBlockSum (Matrix, Matrix) := Matrix => (A,B) -> (
    if numColumns A == 0 then return B;
    if numColumns B == 0 then return A;    
    return A ++ B
    
)

-- Basic functions for manipulating matrices we may want to use

-- Check if a matrix is square
isSquare = method()
isSquare (Matrix) := Boolean => M -> (
    numRows(M) == numColumns(M)
)

-- Check if a matrix is square and symmetric
isSquareAndSymmetric = method()
isSquareAndSymmetric (Matrix) := Boolean => M -> (
    transpose(M) == M
)

-- Check if a matrix represents a degenerate bilinear form
isDegenerate = method()
isDegenerate (Matrix) := Boolean => M ->(
    det(M) == 0
    )

-- Check if a square matrix is upper left triangular
isUpperLeftTriangular = method()
isUpperLeftTriangular (Matrix) := Boolean => M -> (

    if not isSquare(M) then error "Error: matrix isn't square";

    n := numRows(M);
    
    for i from 0 to n-1 do(
	for j from 0 to n-1 do(
	    
	    -- Search in the matrix entries that lie below the main antidiagonal
	    if i + j >= n then(
		
		-- If any entry in this range is nonzero then the matrix isn't upper left triangular
		if not M_(i,j) == 0 then(
		    return false
		    );
		);
        );
    );
true
)
 

-- Check if a square matrix is diagonal
isDiagonal = method()
isDiagonal (Matrix) := Boolean => M -> (

    if not isSquare(M) then error "Error: matrix isn't square";

    n := numRows(M);
    
    for i from 0 to n-2 do(
	for j from  i+1 to n-1 do(
	    
	    -- Search in the matrix entries that aren't on diagonal
	    if i != j  then(
		
		-- If any entry off diagonal  is nonzero then the matrix isn't diagonal
		if  M_(i,j) != 0 or M_(j,i) != 0  then(
		    return false
		    );
		);
        );
    );
true
)


-- TODO we should rename this if we want to export it.
--Diagonalize method
--Given a symmetric matrix, this function outputs a diagonal matrix congruent to original matrix
diagonalize = method()
diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    k := ring AnonMut;
    if isField k == false then error "Error: expected matrix entries from a field";
    A := mutableMatrix AnonMut;
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    n := numRows(A);
    for col from 0 to (n-1) do (
	--If diagonal entry in column "col" is zero
        if A_(col,col) == 0 then (
            for row from col+1 to n-1 do ( 
		--Scan for nonzero entries below the diagonal entry
                if A_(row,col) != 0 then (
                    if A_(row,row) == 0 then (
		        --Row reduction to make A_(col,col) nonzero
                        rowAdd(A,col,1,row);
		        --Column reduction to keep reduced matrix congruent to original matrix
                        columnAdd(A,col,1,row);
                        )
                    else (
		        --Row and column swaps to make A_(col,col) nonzero
                        rowSwap(A,col,row);
                        columnSwap(A,col,row);
                        );
                    break;
                    );
                );
            );
        --Now A_(col,col) != 0 unless there was a zero row/column and we use it to clear the column below
        if A_(col,col) != 0 then (
            for row from (col+1) to (n-1) do (
                temp:=A_(row,col);
                --More row reduction make every entry below A_(col,col) is zero
                rowAdd(A,row,-temp/A_(col,col),col);
	        --Column reduction to keep reduced matrix congruent
                columnAdd(A,row,-temp/A_(col,col),col);
                );
            );
        );
    return matrix A 
    )

--diagonalizeOverInt method
--given a symmetric matrix, this function outputs a diagonal matrix congruent to original matrix
--capable of diagonalizing over rings (algorithm has no divisions)
diagonalizeOverInt = method()
diagonalizeOverInt (Matrix) := (Matrix) => (AnonMut) -> (
    A := mutableMatrix AnonMut;
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    n:=numRows(A);
    for col from 0 to (n-1) do (
        if A_(col,col) == 0 then (
            for row from col+1 to n-1 do ( 
		--Scan for nonzero entries below the diagonal entry
                if A_(row,col) != 0 then (
                    if A_(row,row) == 0 then (
		        --Row reduction to make A_(col,col) nonzero
                        rowAdd(A,col,1,row);
		        --Column reduction to keep reduced matrix congruent to original matrix
                        columnAdd(A,col,1,row);
                        )
                    else (
		        --Row and column swaps to make A_(col,col) nonzero
                        rowSwap(A,col,row);
                        columnSwap(A,col,row);
                        );
                    break;
                    );
                );
            );
        --Now A_(col,col) != 0 unless there was a zero row/column and we use it to clear the column below
        if A_(col,col) != 0 then (
            for row from (col+1) to (n-1) do (
                temp:=A_(row,col);
                rowMult(A,row,A_(col,col)); --multiply row row by A_(col,col)
                columnMult(A,row,A_(col,col)); --column multiplication to keep reduced matrix congruent
                rowAdd(A,row,-temp,col); --more row reduction make every entry below A_(col,col) is zero
                columnAdd(A,row,-temp,col); --column reduction to keep reduced matrix congruent
                );
            );
        );
    return matrix A 
    )

-- This function aims to find the radical of a quadratic space.
-- This is reliant on the diagonalize() method being applicable for singular matrices

truncateRadical=method()
truncateRadical(Matrix):=(Matrix)=> (A) -> (
    truncatedMatrix:= mutableMatrix A;
    if not numRows(A)==numColumns(A) then (
        print ("Input is not a square matrix");
    )
    else (
        truncatedMatrix=mutableMatrix diagonalize(A);
        foundRadical:=false;
        for i from 0 to (numRows(A)-1) do (
            if truncatedMatrix_(i, i)==0 then (
                foundRadical=true;
                break
            );
            print("The quadratic space does not have a radical!");
        );
            if foundRadical===true then (
                n:=numRows(A)-1;
                for i from 0 to n do (
                truncatedMatrix=mutableMatrix submatrix'(matrix truncatedMatrix, {i}, {i});
                    if (n>0) then (n=n-1;)
                    else (break);
            );
            B:=matrix truncatedMatrix;
            return B;
        );
    );
);





--given diagonal matrix, split off any <a>+<-a> and return number of times we can do this as well as smaller matrix with none of these

splitOffObviousHyperbolic = method()
splitOffObviousHyperbolic (Matrix) := (ZZ,Matrix) => (A) -> (
    --matrix must be symmetric
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    foundHyperbolic := 0;
    remainingMatrix := A;
    for i from 0 to (numRows(A)-1) do (
        for j from (i+1) to (numRows(A)-1) do (
            if (A_(i,i) == -A_(j,j) and A_(i,i) != 0) then (
                 foundHyperbolic = 1;
                 remainingMatrix = submatrix'(A,{i,j},{i,j});
                 return(foundHyperbolic,remainingMatrix)
                 );
            );
        );
    return(foundHyperbolic,remainingMatrix);
    )


splitOffObviousHyperbolics = method()
splitOffObviousHyperbolics (Matrix) := (ZZ,Matrix) => (A) -> (
    numberHyperbolics:=0;
    notFinished:=1;
    while (notFinished == 1) do (
        currentState:= splitOffObviousHyperbolic(A);
        notFinished = currentState_0;
        numberHyperbolics = numberHyperbolics + notFinished;
        A = currentState_1;
        );
    return (numberHyperbolics,A);
    )

-- Takes in symmetric matrix over QQ and diagonalizes, removes squares from entries, and splits off hyperbolic forms that immediately appear as <a> + <-a>
rationalSimplify = method()
rationalSimplify (Matrix) := (ZZ,Matrix) => (A) -> (
    --matrix must be symmetric
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    --diagonalize the matrix
    B:= mutableMatrix(diagonalize(A));
    --replace entry with smallest magnitude integer in square class
    for i from 0 to (numRows(B)-1) do (
        B_(i,i) = squarefreePart(B_(i,i));
        );
    C := matrix B;
    --split off hyperbolic forms <a> + <-a>
    return(splitOffObviousHyperbolics(C))
    )



--Nikita Borisov and Frenly Espino

-- Given a matrix A, wittDecomp decomposes A as a sum nH + Q, where n= number of hyperbolic forms, and Q is (presumed to be) anisotropic.  
-- Input:  A square symmetric matrix A over an exact field (not RR or CC)
-- Output: (n, Q), n=integer giving number of hyperbolic spaces, Q=anisotropic forms.

-- Future Work: Should be able to input a bound.

wittDecomp = method()
wittDecomp (Matrix) := (ZZ,Matrix) => (A) -> (
    k:= ring A;   
  
    -- Add error in case the base field is RR or CC
    if (instance(k,InexactFieldFamily) or instance(k,RealField) or instance(k,ComplexField)) then error "Error: base field is inexact, use wittDecompInexact() instead";
    
    n:=numRows(A); --rank of matrix
    x:= symbol x;
    R:=k[x_0..x_(n-1)]; -- local variable ring
    
    -- Create quadratic from f from matrix A
    f:=sum (
        for i from 0 to (n-1) list (
            sum (for j from 0 to (n-1) list (A_(i,j)*x_i*x_j))
            )
        );
 
    -- will use rationalPoints package to see if f has a zero.  If so, A is isotropic.  
    -- rationalPoints package seeks a zero upto variable bound
    
    use k;
    solnPt:=new MutableList;
    solnFound:= false;
    for bound from 1 to 10 do (
        solns:= new MutableList from rationalPoints(ideal(f),Bound=>bound);
        if (#solns >1) then(
            if solns#0 == toList(n:0) then (solnPt=solns#1) else (solnPt=solns#0);
            solnFound =true;
            break;
        );
    );


    --if no solutions found we assume the form is anisotropic
    if ( not solnFound) then (return (0,A));
    
    --if solution found for rank 2 form, then the form is purely hyperbolic
    if ((n==2) and  (det(A)==0))then (error "Matrix singular, run wittDecompGeneral instead" );
    if (n==2) then (return (1,matrix(k,{{}})));

    -- if found a solution, record it as row matrix z. Then find vector y such that <z,y> <>0.  
    -- z and y then generate a hyberbolic plane. 
    z:=matrix{toList solnPt}; --z as a row matrix
    zA:=z*A; --z*A
    y :=new MutableMatrix from matrix{toList(n:(0/1))};
    for i from 0 to (n-1) do (
        if (zA_(0,i) != 0) then (y_(0,i)=1; break;);
    );
    
    --now z and y span a copy of |H in the bilinear form
    --we need to find a basis of vectors orthogonal (wrt bilinear form) to x and y 
    orthoComp :=gens kernel((z||matrix(y))*A);   -- a (n-2) x n matrix.  
    
    --now recursively apply wittDecomp to orthoComp^T*A*orthoComp a (n-2)-by-(n-2) Gram matrix
    subComputation := wittDecomp(transpose(orthoComp)*A*orthoComp);
    
    -- subComputation_0 gives number of hyperbolic forms in (n-2) x (n-2) subform.  
    -- 1+ subComputation_0 is the number of hyperbolic forms in A
    -- subComputation_1 is the anisotropic part of A and (also) the subform part.
    
    return (1+subComputation_0, subComputation_1);
)


--wittDecomp method for InexactFieldFamily

-- wittDecompInexact calculates the Witt decomposition for matrix A over the fields kk=RR, CC
-- Function assumes that A has maximal rank. 

wittDecompInexact=method()
wittDecompInexact (Matrix) := (ZZ,Matrix) => (A) -> (
    k := ring A;
    
    -- checks that we have RR or CC as our field
    if not (instance(k,RealField) or instance(k,ComplexField) or instance(k,InexactFieldFamily)) then error "Error: base field is not RR or CC";
    
    n:=numRows(A); --rank of matrix
    
    --if k is the complex numbers, witt decomposition depends only on rank
    if (k===CC or instance(k,ComplexField)) then (
        if (n%2==0) then(return (n//2,matrix(CC,{{}}))) --if rank is even, then matrix decomposes into n/2 hyberbolic forms with no anisotropic parts
        else return (n//2,id_(k^1)); --if rank is odd, matrix decomposes into (n-1)/2 hyperbolic forms with 1-by-1 anisotropic part
        );
    
    --if k is the real numbers, witt decomposition depends on rank and signature
    if (k===RR or instance(k,RealField)) then (
        diagA := diagonalize(A);
        posEntries := 0; --for loop counts the number of positive diagonal entries of diagA
        negEntries := 0; --for loop counts the number of negative diagonal entries
	for i from 0 to (n-1) do(
            if diagA_(i,i)>0 then(
                posEntries=posEntries+1;
            );
	    if diagA_(i,i)<0 then(
                negEntries=negEntries+1;
            );
        );

        if (posEntries + negEntries > n) then (print "A is singular");
        wittIndex := min(posEntries,negEntries); -- witt index is given by how many positive-negative diagonal entry pairs exist
        signature := posEntries-negEntries; 
        if signature == 0 then (return (wittIndex,matrix(RR,{{}})))
        else if signature > 0 then ( return (wittIndex, id_(k^(signature)))) --signature characterizes anisotropic part
        else return (wittIndex, -id_(k^(-signature)));
        );
);








---------------
-- Types
---------------

-- We define GrothendieckWittClass to be a new type,
--    meant to represent the isomorphism class of a symmetric bilinear form
--    over a base field.
GrothendieckWittClass = new Type of HashTable
GrothendieckWittClass.synonym = "Grothendieck Witt Class"

gwClass = method()

-- A class in GW can be constructed from a representing matrix
gwClass (Matrix) := GrothendieckWittClass => M -> (
  new GrothendieckWittClass from {
      symbol matrix => M,
      symbol cache => new CacheTable
      }
)

-- This allows us to extract the matrix from a class
matrix GrothendieckWittClass := Matrix => beta -> beta.matrix

-- Check if a constructed class is well-defined
isWellDefined GrothendieckWittClass := Boolean => beta -> (
    
    -- Returns false if a matrix isn't symmetric
    --	  Note that this will also return false if a matrix isn't square,
    --	      so we don't need another check for that.
    if not (transpose(beta.matrix) == beta.matrix) then(
	<< "-- Defining matrix is not symmetric" << endl;
	return false
	);
    
    -- Returns false if a matrix is defined over a field of characteristic two
    if char ring beta.matrix == 2 then(
	<< "-- Package does not support base fields of characteristic two" <<endl;
	return false
	);
    
    -- Returns false if the matrix has determinant zero
    if det beta.matrix == 0 then(
	<< "-- Matrix defines a degenerate bilinear form" <<endl;
	return false
        );
    
    -- Returns false if the defining matrix isn't defined over a field
    if not isField ring beta.matrix then(
	<< "-- Matrix is not defined over a field" <<endl;
	return false
	);
    
    true);

-- Method for returning the ring the matrix is defined over
baseField = method()
baseField GrothendieckWittClass := Ring => beta -> (
    if(isWellDefined(beta)===true) then (ring beta.matrix)
)

-- Method for taking a direct sum of two Grothendieck-Witt classes
gwAdd = method()

gwAdd(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (
    
    -- Returns an error if the underlying fields of the two inputted symmetric bilinear forms are different
    if not ring beta.matrix === ring gamma.matrix then error "Error: these classes have different underlying fields";

    b := beta.matrix;
    g := gamma.matrix;
    
    return gwClass(b++g)
    )

-- Method for taking a tensor product of two Grothendieck-Witt classes
gwMultiply = method()

gwMultiply(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (

    -- Returns an error if the underlying fields of the two inputted symmetric bilinear forms are different
    if not ring beta.matrix === ring gamma.matrix then error "Error: these classes have different underlying fields";

    b := beta.matrix;
    g := gamma.matrix;
    
    return gwClass(b**g)
    )




------------
-- Simplifying forms
------------

---------
-- diagonalForm method
-- inputs a GWClass and outputs its diagonal form as a GWClass
---------

diagonalForm = method()
diagonalForm (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    
    -- TODO: quick check if the form is already diagonal, if so return it and cache it
    
    
    -- Check if the diagonalForm has already been computed, if so recall it from the cache
    if beta.cache.?diagonalForm then return beta.cache.diagonalForm;
    
    -- Number of rows of beta
    n := numRows(beta.matrix);
    
    -- If the field is the complex numbers, its diagonal form is an identity matrix of the same size
    if (baseField(beta) === CC or instance(baseField(beta),ComplexField)) then(
	identityMat := matrix(mutableIdentity(CC,n));
	
	-- Cache the answer before returning
	beta.cache.diagonalForm = gwClass(identityMat);
	return gwClass(identityMat)
	);
    
    -- If the field is the real numbers, we can run wittDecompInexact to determine its Witt index and anisotropic part
    if (baseField(beta) === RR or instance(baseField(beta),RealField)) then(
	
	-- Get wittIndex and anisotropic part
	(wittIndex,anisotropicPart) := wittDecompInexact(beta.matrix);
	
	-- Make empty output matrix to populate
	diagOutputMatrix := matrix(RR,{{}});
	H := matrix(RR,{{1,0},{0,-1}});
	
	
	
	-- Sum as many hyperbolic forms as the Witt index
	
	if wittIndex > 0 then(
	    for i in 1..wittIndex do(
		diagOutputMatrix = safeBlockSum(diagOutputMatrix,H);
		);
	    
	    );
	
	-- Add on the anisotropic part
	diagOutputMatrix = safeBlockSum(diagOutputMatrix, anisotropicPart);
	
	-- Cache and return
	beta.cache.diagonalForm = diagOutputMatrix;
	return gwClass(diagOutputMatrix)
	
	);
    
    betaMatrix := beta.matrix;
    diagonalFormOfBetaMatrix := diagonalize(betaMatrix);
    beta.cache.diagonalForm = gwClass(diagonalFormOfBetaMatrix);
    return gwClass(diagonalFormOfBetaMatrix) 
    );




-- simplifyForm method

simplifyFormVerbose = method()
simplifyFormVerbose (GrothendieckWittClass) := (GrothendieckWittClass, String) => beta -> (
    -- Check if the diagonalForm has already been computed, if so recall it from the cache    
    gamma := diagonalForm(beta);
    
    -- Get the matrix, base field, and rank of the form
    A := gamma.matrix;
    k := baseField(gamma);
    n := numRows(gamma.matrix);
    
    -- If the form has no rank, stop here and return the zero form. This gets rid of fringe cases later on
    if n == 0 then return (gamma,"zero form");
    
    -- Define the hyperbolic form over our field
    H := matrix(k,{{1,0},{0,-1}});
        
    -----------------------
    -- Field by field cases
    -----------------------
        
    --In the case of C:
    if (k === CC or instance(k,ComplexField)) then(
        if even n then(
        
            -- The form is (n/2)H
            numHypFormsCCEven := n/2;
        
            simplifiedFormCCEven := H;
        
            for i in 1..sub((numHypFormsCCEven-1),ZZ) do(
                simplifiedFormCCEven = simplifiedFormCCEven ++H;
		);
        
            outputStringCCEven := toString(numHypFormsCCEven) | "H";
            return (gwClass(simplifiedFormCCEven),outputStringCCEven)
            );
        if odd n then(
        
            -- The form is floor(n/2)H + <1>
            numHypFormsCCOdd := floor(n/2);
            simplifiedFormCCOdd := H;
        
            for i in 1..sub((numHypFormsCCOdd-1),ZZ) do(
                simplifiedFormCCOdd = simplifiedFormCCOdd ++H;
		);
        
            simplifiedFormCCOdd = simplifiedFormCCOdd ++ matrix(k,{{1}});
        
            outputStringCCOdd := toString(numHypFormsCCOdd) | "H + <1>";
            return (gwClass(simplifiedFormCCOdd),outputStringCCOdd)
            );
	);
	-- End CC case
    
    -- If the field is R, look at the sign of elements along the diagonal
    if (k === RR or instance(k,RealField)) then(
	print("field is R");
        posEntries := 0; --for loop counts the number of positive diagonal entries of diagA
        negEntries := 0; --for loop counts the number of negative diagonal entries
	for i from 0 to (n-1) do(
            if A_(i,i)>0 then(
                posEntries=posEntries+1;
            	);
            if A_(i,i)<0 then(
                negEntries=negEntries+1;
            	);
            );
            
	-- Number of hyperbolic forms is the number of positive and negative entries
        wittIndexRR := min(posEntries,negEntries);
    	
	-- Make an empty matrix and string to add output to
	simplifiedFormRR := matrix(k,{{}});
	outputStringRR := "";
	
	-- Add hyperbolic forms to output
	for i in 1..(wittIndexRR) do(
	    simplifiedFormRR = safeBlockSum(simplifiedFormRR,H);
            );
	
	if wittIndexRR> 0 then(
	    outputStringRR = outputStringRR | toString(wittIndexRR) | "H";
	    );
	
	-- Look at what's left over
	posEntries = posEntries - wittIndexRR;
	negEntries = negEntries - wittIndexRR;
	
	-- Build the anisotropic part
	anisotropicPart := safeBlockSum(matrix(mutableIdentity(k,posEntries)),((-1)*matrix(mutableIdentity(k,negEntries))));
	
	-- Add on (safely) the anisotropic part
	simplifiedFormRR = safeBlockSum(simplifiedFormRR,anisotropicPart);
	
	-- Amend the string output with anisotropic info
	if posEntries > 0 then(
	    outputStringRR = outputStringRR | " + " | toString(posEntries) | "<1>";
	    );
	
	if negEntries > 0 then(
	    outputStringRR = outputStringRR | " + " | toString(negEntries) | "<-1>";
	    );
	
	return (gwClass(simplifiedFormRR), outputStringRR);
	
	);
	-- End RR case
	-- TODO are we getting the string "4H" or whatever in the output string?
       
       
    -- Finite field case
    
    if instance(k, GaloisField) then(
	
	-- Take the diagonal entries to a list
	diagEntries := {};
	for i from 0 to sub((n-1),ZZ) do(
	    diagEntries = append(diagEntries, A_(i,i));
	    );
	
	-- Start counting squares and nonsquares
        numSquares:= 0;
	numNonSquares:= 0;
	
	-- Define a new variable which we will set equal to a nonsquare representative
	nonSquareRepresentative := sub(0,k);
	
	
	-- Check if each elt is a square or not
	for x in diagEntries do(
	    
	    -- If it is a square
	    if legendreBoolean(x) then(
		numSquares = numSquares + 1;
		);
	    
	    -- If it is not a square
	    if not legendreBoolean(x) then(
		--Pick x to be a representative of the non-square class in F_q^x / (F_q^x)^2
		nonSquareRepresentative = x;
		numNonSquares = numNonSquares + 1;
		);
	    );
	
	
	-- If -1 is a square class in k, then 2<1> = H
	if legendreBoolean(sub(-1,k)) then(
	    
	    -- Set up output matrix and string
	    simplifiedFormGFSquare := matrix(k,{{}});
	    outputStringGFSquare := "";	    
	    
	    -- Number of hyperbolic forms
            wittIndexGFSquare := floor(numSquares/2) + floor(numNonSquares/2);
	    
	    -- If there are hyperbolic forms add them to the matrix and string
	    if wittIndexGFSquare > 0 then(
		outputStringGFSquare = outputStringGFSquare | toString(wittIndexGFSquare) | "H";
		for i in 1..(wittIndexGFSquare) do(
		    simplifiedFormGFSquare = safeBlockSum(simplifiedFormGFSquare,H);
		    );
		);
	    
	    -- There will be at most one extra square or one extra nonsquare left over	
	    if odd numSquares then(
		outputStringGFSquare = outputStringGFSquare | " + <1>";		
		simplifiedFormGFSquare = safeBlockSum(simplifiedFormGFSquare ,matrix(k,{{1}}));		
		);
	    
	    if odd numNonSquares then(
		outputStringGFSquare = outputStringGFSquare | " + <" | toString(nonSquareRepresentative) | ">";		
		simplifiedFormGFSquare = safeBlockSum(simplifiedFormGFSquare,((-1)*matrix(k,{{nonSquareRepresentative}})));
		
		);
		
	    return(gwClass(simplifiedFormGFSquare),outputStringGFSquare)
	    
	    
	    ); -- End "-1 is a square in GF"
	
	-- If -1 is not a square then it is a representative for the nonsquare class	  	   
	   if not legendreBoolean(sub(-1,k)) then (
	      	            
	       -- Set up output matrix and string
	       simplifiedFormGFNonSquare := matrix(k,{{}});
	       outputStringGFNonSquare := "";
	       
	       
	       -- Number of hyperbolic forms
	       wittIndexGFNonSquare := min(numSquares,numNonSquares);	       
	       
	       -- Look at any remaining squares or nonsquares. One of these must be zero
	       numSquares = numSquares - wittIndexGFNonSquare;
	       numNonSquares = numNonSquares - wittIndexGFNonSquare;
	       
	       
	       -- Add on hyperbolic part
	       if wittIndexGFNonSquare > 0 then(
		   for i in 1..(wittIndexGFNonSquare) do(
		       simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, H);
		       
		       );
		   outputStringGFNonSquare = outputStringGFNonSquare | toString(wittIndexGFNonSquare) | "H";
		   
		   );
	       
	       -- Add any anisotropic part
	       if numSquares > 0 then(
		   simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, matrix(mutableIdentity(k,numSquares)));
		   outputStringGFNonSquare = outputStringGFNonSquare | " + " | toString(numSquares) | "<1>";
		   );
		   
	      if numNonSquares > 0 then(
		   simplifiedFormGFNonSquare = safeBlockSum(simplifiedFormGFNonSquare, ((-1)*matrix(mutableIdentity(k,numNonSquares))));
		   outputStringGFNonSquare = outputStringGFNonSquare | " + " | toString(numNonSquares) | "<-1>";
		   );
	       
	       return (gwClass(simplifiedFormGFNonSquare), outputStringGFNonSquare)		   
	       
	       ); -- End "-1 is not a square in GF"
	
	    
	  
			
      
	   
	   
	   );-- End finite field case
       
          -- Number field case
    
    -- NOTE: I wanted to do member(QQ,k.baseRings) but this returns false.
    -- 	   This is because baseField(..) will return a ring
    
    
    -- We're making the number field case a catch-all case now
    
    
--    if member(QQ, k.baseRings) and (isField(k)) then(
	
	
	-- Set up output form and matrix
	simplifiedFormQQ := matrix(k,{{}});
        outputStringQQ := "";
	
	-- Get number of confirmed hyperbolic forms and remainder from wittDecomp
	(numHypForms,B) := wittDecomp(A);
	
	-- Add any hyperbolic forms if they exist
	if numHypForms > 0 then(
	    for i in 1..(numHypForms) do(
		
		simplifiedFormQQ = safeBlockSum(simplifiedFormQQ, H);
		outputStringQQ = outputStringQQ | toString(numHypForms) | "H + potentially anisotropic part";
		
		);
	    );
	
	if numHypForms == 0 then(
	    outputStringQQ = outputStringQQ | " potentially anisotropic part";
	    
	    );
	 simplifiedFormQQ = safeBlockSum(simplifiedFormQQ, B);
	 
	 return (gwClass(simplifiedFormQQ),outputStringQQ)
		
       
	
--	);-- End number field case
    
    
    
    );


simplifyForm = method()
simplifyForm (GrothendieckWittClass) := (GrothendieckWittClass) => beta -> (
    beta.cache.diagonalForm = (simplifyFormVerbose(beta))_0;
    return (simplifyFormVerbose(beta))_0
);

simplifyFormString = method()
simplifyFormString (GrothendieckWittClass) := (String) => beta -> (
    return (simplifyFormVerbose(beta))_1
);


----------------------------
-- A1-Brouwer degree methods
----------------------------

globalA1Degree = method()
globalA1Degree (List) := (GrothendieckWittClass) => (Endo) -> (
    -- Endo is the list {f_1, f_2, ..., f_n} of polynomials kk^n -> kk^n
    
    -- n is the number of polynomials
    n := #Endo; -- n is the no. of polynomials    

    -- Get the underlying field, assert it is a field
    kk := coefficientRing(ring(Endo#0)); 
    
    if isField(kk) == false then(
	kk = toField(kk);
	);
    
    -- Let S = k[x_1..x_n] be the ambient polynomial ring
    S:=ring(Endo#0);
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then error "Error: morphism does not have isolated zeroes";
    
    -- Check the number of variables matches the number of polynomials
    if not #(gens S) == n then error "Error: the number of variables does not match the number of polynomials.";
    
    -- If the field is CC, just output gwClass of an identity matrix of rank = rankAlgebra
    if (kk === CC or instance(kk,ComplexField)) then(
    	rankAlgebra:=rankGlobalAlgebra(Endo);
    	return gwClass(matrix(mutableIdentity(CC,rankAlgebra)));
    );
    
    -- If the field is RR, ask the user to run it over QQ instead, then simplify over RR
    if (kk === RR or instance(kk,RealField)) then error "Error: globalA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output a GrothendieckWittClass. Then extract the matrix, base change it to RR, and then run simplifyForm().";    
    
    -- Create internal rings/matrices
    
    -- Initialize a polynomial ring in X_i's and Y_i's to compute the Bezoutian in
    X:=local X;
    Y:=local Y;
 --   Rtemp := kk(monoid[X_1..X_n]);
 --   R := Rtemp(monoid[Y_1..Y_n]);
 
   R:=kk(monoid[X_1..X_n|Y_1..Y_n]);
    -- Create an (n x n) matrix D which will be populated by \Delta_{ij} in the paper
    D := "";
    try D = mutableMatrix id_((frac R)^n) else D= mutableMatrix id_(R^n);
    
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(
	    -- iterate through the entries of the matrix D and populate it with the following information ...
        -- create the list {Y_1, ..., Y_(j-1), X_j, ..., X_n}. Note Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList1 := apply(toList (Y_1..Y_j|X_(j+1)..X_n),i->i_R);
        -- create the list {Y_1, ..., Y_j, X_(j+1), ..., X_n }. Note Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList2 := apply(toList (Y_1..Y_(j+1)| X_(j+2)..X_n),i->i_R); 
        -- suppose our endomorphisms are given in the variables x_1, ..., x_n.
        -- map f_i(x_1, ..., x_n) to f_i(Y_1, ..., Y_(j-1), X_j, ..., X_n) resp.
        -- then take the difference f_i(Y_1, ..., Y_(j-1), X_j, ..., X_n) - f_i(Y_1, ... Y_j, X_(j+1), ..., X_n)
        numeratorD := ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i)); 
        -- divide it by X_j - Y_j, Note Macaulay2 is 0-indexed hence the difference in notation. 
	    D_(i,j)= numeratorD/((X_(j+1))_R-(Y_(j+1))_R); 
	); 
    );
    
    -- Set up the local variables bezDet and bezDetR
    bezDet:="";
    bezDetR:="";


    -- The determinant of D is interpreted as living in Frac(k[x_1..x_n]),
    -- so we can try to lift it to k[x_1..x_n]           
    if liftable(det(D),R) == true then(
	bezDetR = lift(det(D),R);
	);
    
    -- In some computations, applying lift(-,R) won't work. So we instead lift
    -- the numerator and then divide by a lift of the denominator (which will be
    -- a scalar) to the coefficient ring k
    if not liftable(det(D),R) == true then(
	bezDet = lift(numerator(det(D)), R) / lift(denominator(det(D)),coefficientRing R);
    	bezDetR = lift(bezDet, R);
	);

    -- Define formal variables X_i, Y_i that replace x_i
    RX:=kk[X_1..X_n]; 
    RY:=kk[Y_1..Y_n];

    -- mapxtoX replaces all instances of x_i with X_i. mapxtoY does the same but with Y_i's
    mapxtoX:= (map(RX,S,toList(X_1..X_n))); 
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n)));

    -- Compute the standard basis of kk[X_1, ..., X_n]/(f_1, ..., f_n)
    standBasisX := basis (RX/(ideal (leadTerm (mapxtoX ideal Endo)))); 
    standBasisY := basis (RY/(ideal (leadTerm (mapxtoY ideal Endo)))); 

    -- defines an ideal (f_1(X), ..., f_n(X))
    id1 := (ideal apply(toList(0..n-1), i-> mapxtoX(Endo_i))); 
    -- defines an ideal (f_1(Y), ..., f_n(Y))
    id2 := (ideal apply(toList(0..n-1), i-> mapxtoY(Endo_i))); 

    -- takes the sum of the ideals (f_1(X),...,f_n(X)) + (f_1(Y),...,f_n(Y)) in the ring kk[X_1..Y_n]
    promotedEndo :=sub(id1,R)+sub(id2,R); 

    -- Here we're using that (R/I) \otimes_R (R/J) = R/(I+J) in order to express Q(f) \otimes Q(f), where X's are the variables in first term, Y's are variables in second par
    Rquot:= R/promotedEndo; 


    -- moves the standard bases to the quotient ring
    sBXProm := sub(standBasisX, Rquot); 
    sBYProm := sub(standBasisY, Rquot); 
    
    --reduces the bezDetR determinant subject to the ideal in the X's and Y's
    bezDetRed := bezDetR % promotedEndo;

    
    -- define a ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring
    phi0 := map(kk,Rquot,(toList ((2*n):0))); 

    -- m is the dimension of the basis for the algebra
    m:= numColumns(sBXProm);

    -- Now create Bezoutian matrix B for the quadratic form by reading off the coefficients. 
    -- B is a (m x m) matrix.  Coefficent B_(i,j) is the coefficient of the (ith basis vector x jth basis vector) in tensor product.
    -- phi0 maps the coefficient to kk
    B:= mutableMatrix id_(kk^m);
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            B_(i,j)=phi0(coefficient((sBXProm_(0,i)**sBYProm_(0,j))_(0,0), bezDetRed));
        );
    );
    return gwClass(matrix(B));

);


-----
-- Local
-----
localA1Degree = method()
localA1Degree (List, Ideal) := (GrothendieckWittClass) => (Endo,p) -> (
    -- Endo is the list {f_1, f_2, ..., f_n} of polynomials kk^n -> kk^n
    
    -- n is the number of polynomials
    n := #Endo; -- n is the no. of polynomials
    

    -- Get the underlying field, assert it is a field
    kk := coefficientRing(ring(Endo#0)); 
    
       
    if isField(kk) == false then(
	kk = toField(kk);
	);

    
    -- Let S = k[x_1..x_n] be the ambient polynomial ring
    S:=ring(Endo#0);
    
    -- Create the ideal J.  It has prop that  k[x_1 .. x_n]_p/I_p is isom to k[x_1 .. x_n]/J    
    J:=(ideal Endo):saturate(ideal Endo,p);
    
    -- Get the dimension of the local algebra Q_p(f)
    localFormRank := numColumns(basis(S/J));
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then error "Error: morphism does not have isolated zeroes";
    
    -- Check the number of variables matches the number of polynomials
    if not #(gens S) == n then error "Error: the number of variables does not match the number of polynomials.";
    
    -- If the field is CC, just output gwClass of an identity matrix of rank = localFormRank
    if (kk === CC or instance(kk,ComplexField)) then(
	return gwClass(matrix(mutableIdentity(CC,localFormRank)))
    );
    
    -- If the field is RR, ask the user to run it over QQ instead, then simplify over RR
    if (kk === RR or instance(kk,RealField)) then error "Error: localA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output a GrothendieckWittClass. Then extract the matrix, base change it to RR, and then run simplifyForm().";
    
    
    -- Create internal rings/matrices
    
    -- Initialize a polynomial ring in X_i's and Y_i's to compute the Bezoutian in
    X := local X;
    Y := local Y;
   -- R := kk[X_1..Y_n];
    R := kk(monoid[X_1..X_n|Y_1..Y_n]);
    
    -- Create a matrix D which will be populated by \Delta_{ij} in the paper
    D := "";
    try D = mutableMatrix id_((frac R)^n) else D= mutableMatrix id_(R^n);
    
    for i from 0 to (n-1) do (
	for j from 0 to (n-1) do(
	    -- iterate through the entries of the matrix D and populate it with the following information ...
        -- create the list {Y_1, ..., Y_(j-1), X_j, ..., X_n}. Note Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList1 := apply(toList (Y_1..Y_j|X_(j+1)..X_n),i->i_R);
        -- create the list {Y_1, ..., Y_j, X_(j+1), ..., X_n }. Note Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList2 := apply(toList (Y_1..Y_(j+1)| X_(j+2)..X_n),i->i_R); 
        -- suppose our endomorphisms are given in the variables x_1, ..., x_n.
        -- map f_i(x_1, ..., x_n) to f_i(Y_1, ..., Y_(j-1), X_j, ..., X_n) resp.
        -- then take the difference f_i(Y_1, ..., Y_(j-1), X_j, ..., X_n) - f_i(Y_1, ... Y_j, X_(j+1), ..., X_n)
        numeratorD := ((map(R,S,targetList1))(Endo_i)-(map(R,S,targetList2))(Endo_i)); 
        -- divide it by X_j - Y_j, Note Macaulay2 is 0-indexed hence the difference in notation. 
	    D_(i,j)= numeratorD/((X_(j+1))_R-(Y_(j+1))_R); 
	); 
    );
    
    -- Set up the local variables bezDet and bezDetR
    bezDet:="";
    bezDetR:="";
    
    
    -- The determinant of D is interpreted as living in Frac(k[x_1..x_n]),
    -- so we can try to lift it to k[x_1..x_n]           
    if liftable(det(D),R) == true then(
	bezDetR = lift(det(D),R);
	);
    
    -- In some computations, applying lift(-,R) won't work. So we instead lift
    -- the numerator and then divide by a lift of the denominator (which will be
    -- a scalar) to the coefficient ring k
    if not liftable(det(D),R) == true then(
	bezDet = lift(numerator(det(D)), R) / lift(denominator(det(D)),coefficientRing R);
    	bezDetR = lift(bezDet, R);
	);    
    
    -- Define formal variables X_i, Y_i that replace x_i
    RX:=kk[X_1..X_n]; 
    RY:=kk[Y_1..Y_n];
    
    -- mapxtoX replaces all instances of x_i with X_i. mapxtoY does the same but with Y_i's
    mapxtoX:= (map(RX,S,toList(X_1..X_n)));
    mapxtoY:=(map(RY,S,toList(Y_1..Y_n)));

    -- Find standard basis and define local quotient ring
    list1 := (apply(toList(0..n-1), i-> mapxtoX(Endo_i))); -- list (f_1(X), ..., f_n(X))
    list2 := (apply(toList(0..n-1), i-> mapxtoY(Endo_i))); -- list (f_1(Y), ..., f_n(Y))
    
    -- Apply localAlgebraBasis to compute standard basis for localization
    standBasisX := localAlgebraBasis(list1,mapxtoX p);
    standBasisY := localAlgebraBasis(list2,mapxtoY p); 
    

    
    -- Take the ideal sum of J in the X_i's with J in the Y_i's
    localIdeal :=sub(mapxtoX(J),R)+sub(mapxtoY(J),R);
    
    
    
    
    -- Here we're using that (R/I) \otimes_R (R/J) = R/(I+J) in order to express Q_p(f) \otimes Q_p(f), where X's are the variables in first term, Y's are variables in second par
    Rquot:=R/localIdeal;

    
    
    -- Move the standard bases to the quotient ring
    sBXProm :=apply(toList(0..#standBasisX-1),i-> sub(standBasisX_i,Rquot));
    sBYProm :=apply(toList(0..#standBasisY-1),i-> sub(standBasisY_i,Rquot));
   
    -- Now reduce the bezDetR determinant subject to the local ideal in both the X's, Y's.
    bezDetRed := bezDetR % localIdeal;
    
    -- ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring 
    phi0 := map(kk,Rquot,(toList ((2*n):0))); 
    
    -- m is the dimension of the basis for the local ring
    m:= #sBXProm;
    
    -- Now create Bezoutian matrix B for the quadratic form by reading off the local coefficients. 
    -- B is a (m x m) matrix.  Coefficent B_(i,j) is the coefficient of the (ith basis vector x jth basis vector) in tensor product.
    -- phi0 maps the coefficient to kk
    
    B:= mutableMatrix id_(kk^m);
    for i from 0 to m-1 do (
        for j from 0 to m-1 do (
            B_(i,j)=phi0(coefficient((sBXProm_i**sBYProm_j)_(0,0), bezDetRed));
        );
    );
    return gwClass(matrix(B));
);


---------------------------------------
-- Classifying and comparing GW classes
---------------------------------------


numPosEntries = method()
numPosEntries (GrothendieckWittClass) := ZZ => beta ->(
    B := beta.matrix;
    n := numRows(B);
    kk := ring B;
    if not (kk === RR or instance(kk,RealField) or kk === QQ) then(
        error "Field is not QQ or RR";
        );
    diagB := diagonalize(B);
    posEntries := 0;
    for i from 0 to (numRows(B)-1) do (
        if diagB_(i,i) > 0 then(
            posEntries = posEntries+1;
            );
	);

    return posEntries
);


numNegEntries = method()
numNegEntries (GrothendieckWittClass) := ZZ => beta ->(
    B := beta.matrix;
    n := numRows(B);
    kk := ring B;
    if not (kk === RR or instance(kk,RealField) or kk === QQ) then(
        error "Field is not QQ or RR";
        );
    diagB := diagonalize(B);
    negEntries := 0;
    for i from 0 to (numRows(B)-1) do (
        if diagB_(i,i) < 0 then(
            negEntries = negEntries+1;
            );
	);

    return negEntries
);


signature = method()
signature (GrothendieckWittClass) := ZZ => (beta) ->(
    -- B := beta.matrix;
    -- n := numRows(B);
    -- kk := ring B;
    -- if not (kk === RR or instance(kk,RealField) or kk === QQ) then(
    --     error "Field is not QQ or RR";
    --     );
    -- diagB := diagonalize(B);
    -- posEntries := 0;
    -- negEntries := 0;
    -- for i from 0 to (numRows(B)-1) do (
    --     if diagB_(i,i) > 0 then(
    --         posEntries = posEntries+1;
    --         );
    --     if diagB_(i,i) < 0 then(
    --         negEntries = negEntries+1;
    --         );
    -- 	);
    -- sig := posEntries - negEntries;
    -- return sig
    sig := numPosEntries(beta) - numNegEntries(beta);
    return sig
    );




----
-- Comparing over QQ
----


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
    a:=0;
    if H#?p then (
    	a=H#p;)
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
	    	p := k_i;
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
  
isIsomorphic2 = method()

isIsomorphic2 (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    k1:=baseField(alpha);
    k2:=baseField(beta);
    -- Ensure both base fields are supported
    if not (k1 === CC or instance(k1,ComplexField) or k1 === RR or instance(k1,RealField) or k1 === QQ or instance(k1, GaloisField)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields";
        );
    if not (k2 === CC or instance(k2,ComplexField) or k2 === RR or instance(k2,RealField) or k2 === QQ or instance(k2, GaloisField)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields";
        );
    A:=alpha.matrix;
    B:=beta.matrix;
    -- Ensure both underlying matrices are symmetric
    if A != transpose(A) then (
        error "Underlying matrix is not symmetric";
	);
    if B != transpose(B) then (
        error "Underlying matrix is not symmetric";
	);
    diagA := diagonalize(A);
    diagB := diagonalize(B);
    -- Over CC, diagonal forms over spaces of the same dimension are equivalent if and only if they have the same number of nonzero entries
    if (k1 === CC or instance(k1,ComplexField)) and (k2 === CC or instance(k2,ComplexField)) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        nonzeroEntriesA := 0;
        nonzeroEntriesB := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) != 0 then (
                nonzeroEntriesA = nonzeroEntriesA + 1;
                );
            if diagB_(i,i) != 0 then (
                nonzeroEntriesB = nonzeroEntriesB + 1;
                );
            );
        return (nonzeroEntriesA == nonzeroEntriesB);
        )
    --Over RR, diagonal forms are equivalent if and only if they have the same number of positive, negative, and zero entries
    else if ((k1 === RR or instance(k1,RealField)) and (k2 === RR or instance(k2,RealField))) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        return (signature(alpha)==signature(beta));
        )
    -- Over QQ, call isIsomorphicFormQ, which checks equivalence over all completions
    else if ((k1 === QQ) and (k2 === QQ)) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        return isIsomorphicFormQ(diagA,diagB);
        )
    -- Over a finite field, diagonal forms over spaces of the same dimension are equivalent if and only if they have the same number of nonzero entries and the product of these nonzero entries is in the same square class
    else if (instance(k1, GaloisField) and instance(k2, GaloisField) and k1.order == k2.order) then (
        if (numRows(A) != numRows(B)) then (
            return false;
            );
        countNonzeroDiagA := 0;
        countNonzeroDiagB := 0;
        prodNonzeroDiagA := 1;
        prodNonzeroDiagB := 1;
        for i from 0 to (numRows(A)-1) do (
	    if diagA_(i,i) != 0 then (
		countNonzeroDiagA = countNonzeroDiagA + 1;
                prodNonzeroDiagA = prodNonzeroDiagA * diagA_(i,i);
		);
	    if diagB_(i,i) != 0 then (
		countNonzeroDiagB = countNonzeroDiagB + 1;
                prodNonzeroDiagB = prodNonzeroDiagB * diagB_(i,i);
		);
	    );
        return ((countNonzeroDiagA == countNonzeroDiagB) and (legendreBoolean(prodNonzeroDiagA) == legendreBoolean(prodNonzeroDiagB)));
        )
    -- If we get here, the base fields are not isomorphic
    else error "Base fields are not isomorphic"
    )

--------------------------
-- Checking isotropy
--------------------------

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
	    if (hilbertSymbol(-1, -d, p) == hasseWittInvariant(f, p)) then (
		return false;
		)
	    else (return true);
	    )
	else (
	 -- now use the criteria for n=4
	-- need to check (d=1 and (-1,-1) != =hasseWittInvariant for f) for form to be anisotropic 
	    if (n==4) then (
		if ((equalUptoPadicSquare( d, 1, p)) and  (not (hilbertSymbol(-1,-1,p) == hasseWittInvariant(f, p)))) then (
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
    if (isDegenerate(A)) then (return false);
    
     -- if rank =1, then a non-degenerate form is anisotropic
    if (n==1) then (return true);
    
    --if rank>=5, we can use signature do decide this
    --the non-degenerate form will be anisotropic iff all diagonal entries have same sign
    if (n>= 5) then (
        return ((numPosEntries(alpha) == n) or (numNegEntries(alpha) == n));
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
    disc:= discForm(diagIntEntriesA);
    
    
   
   
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
	    	p := k_i;
		if (isAnisotropicDiagQp(diagIntEntriesA, p)) then (return true);
		i=i+1;
		);

-- if the function hasn't returned false yet, then isotropic over all primes p.  hence form is isotropic over Q	    
  
    return false;
);
	
	
isAnisotropic = method()

isAnisotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    k:=baseField(alpha);
    -- Ensure base field is supported
    if not (k === CC or instance(k,ComplexField) or k === RR or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and k.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    A:=alpha.matrix;
    -- Ensure underlying matrix is symmetric
    if A != transpose(A) then (
        error "Underlying matrix is not symmetric";
	);
    diagA := diagonalize(A);
    -- Over CC, a diagonal form is anisotropic if and only if it is nondegenerate and has dimension 0 or 1
    if (k === CC or instance(k,ComplexField)) then (
        nonzeroEntriesA := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) != 0 then (
                nonzeroEntriesA = nonzeroEntriesA + 1;
                );
            );
        return (nonzeroEntriesA == numRows(A) and numRows(A) <= 1);
        )
    --Over RR, a diagonal form is anisotropic if and only if all of its diagonal entries are positive or all of its diagonal entries are negative
    else if (k === RR or instance(k,RealField)) then (
        posEntriesA := 0;
        negEntriesA := 0;
        for i from 0 to (numRows(A)-1) do (
            if diagA_(i,i) > 0 then (
                posEntriesA = posEntriesA + 1;
                );
            if diagA_(i,i) < 0 then (
                negEntriesA = negEntriesA + 1;
                );
            );
        return ((posEntriesA == numRows(A)) or (negEntriesA == numRows(A)));
        )
    -- Over QQ, call isAnisotropicQ
    else if (k === QQ) then (
        return isAnisotropicQ(alpha);
        )
    -- Over a finite field, a diagonal form is anisotropic if and only if it is nondegenerate, of dimension at most 2, and not the hyperbolic form 
    else if (instance(k, GaloisField) and k.char != 2) then (
        countNonzeroDiagA := 0;
        prodNonzeroDiagA := 1;
        for i from 0 to (numRows(A)-1) do (
	    if diagA_(i,i) != 0 then (
		countNonzeroDiagA = countNonzeroDiagA + 1;
                prodNonzeroDiagA = prodNonzeroDiagA * diagA_(i,i);
		);
	    );
        return ((countNonzeroDiagA == numRows(A)) and (numRows(A) <= 1 or (numRows(A) == 2 and  legendreBoolean(prodNonzeroDiagA) != legendreBoolean(sub(-1,k)))));
        )
    -- We should never get here
    else error "Problem with base field"
    )


isIsotropic = method()
isIsotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    return (not isAnisotropic(alpha));
    )
	    







----------------------------
----------------------------
----------------------------

-- Below here are methods that are NOT USED in local or global a1 degree computations

-- As such I think we should pull them from the final version of the package

----------------------------
----------------------------
----------------------------




----------------------------

-- DOCUMENTATION

beginDocumentation()




document{
    Key => A1BrouwerDegrees,
    Headline => "A package for running A1-Brouwer degree computations in Macaulay2",
    }


undocumented {
    }


document {
    Key => {(squarefreePart, QQ), (squarefreePart, ZZ), squarefreePart},
	Headline => "smallest magnitude representative of a square class over the rationals or integers",
	Usage => "squarefreePart(q)",
	Inputs => {
		QQ => "q" => {"a rational number"},
		-- ZZ => "n" => {"an integer"}
		},
	Outputs => { ZZ => { "the smallest magnitude integer in the square class of ", TT "n"}},
	PARA {"Given a rational number (or integer), ", TT "q", ", this command outputs the smallest magnitude integer, ",
                TEX///$m$///, ", such that ", TEX///$q=lm$///, " for some rational number (or integer) ",
				TEX///$l$///, "."},
	EXAMPLE lines ///
		 squarefreePart(15/72)
		 squarefreePart(-1/3)
	 	 ///,
        }

document{
    Key => GrothendieckWittClass,
    Headline => "a new type, intended to capture the isomorphism class of an element of the Grothendieck-Witt ring of a base field",
    PARA {"A GrothendieckWittClass object is a type of HashTable encoding the isomorphism class of a non-degenerate symmetric bilinear form ", TEX///$V \times V \to k$///, " over a field ", TEX///$k$///, "."},
    PARA{"A GrothendieckWittClass object can be built from a symmetric ", TT "matrix", " over a field using the method ", TT "gwClass()"},
    SeeAlso => {"gwClass"}
    }

document {
    Key => {(gwClass, Matrix), gwClass},
	Headline => "Grothendieck Witt class of a symmetric matrix",
	Usage => "gwClass(M)",
	Inputs => {
		Matrix => "M" => {"a symmetric matrix defined over an arbitrary field"}
		},
	Outputs => { GrothendieckWittClass => { "the isomorphism class of a symmetric bilinear form represented by ", TEX/// $\mathbb{M}$/// }},
	PARA {"Given a symmetric matrix, ", TT "M", ", this command outputs an object of type ", TT "GrothendieckWittClass", ". ",
                "This output has the representing matrix, ", TT "M", ", and the base field of the matrix stored in its CacheTable."},
	EXAMPLE lines ///
		 M := matrix(QQ,{{0,0,1},{0,1,0},{1,0,0}});
		 beta = gwClass(M)
	 	 ///,
	PARA{"The matrix representing a ", TT "GrothendieckWittClass", " element can be recovered using the ", TT "matrix", " command:"},
	EXAMPLE lines ///
	    	beta.matrix
		///,
        PARA{"The base field which the form ", TT "beta", " is implicitly defined over can be recovered with the ", TT "baseField", " method."},
	EXAMPLE lines ///
	    	baseField beta
		///
        }



document {
	Key => {(diagonalize, Matrix), diagonalize},
	Headline => "diagonalizing a symmetric matrix via congruence",
	Usage => "diagonalize(M)",
	Inputs => {
		Matrix => "M" => {"a symmetric matrix over any field"}
		},
	Outputs => { Matrix => { "a diagonal matrix congruent to", TT "M" }},
	PARA {"Given a symmetric matrix ", TT "M", " over any field, this command gives a diagonal matrix congruent to ", TT "M",".",
	"Note that the order in which the diagonal terms appear is not specified."},
	EXAMPLE lines ///
		 M=matrix(GF(17), {{7, 9}, {9, 6}});
		 diagonalize(M)
	 	 ///,
	SeeAlso => {"diagonalForm"}
     	}




document {
    Key => {(baseField, GrothendieckWittClass), baseField},
	Headline => "base field of a Grothendieck Witt class",
	Usage => "baseField(beta)",
	Inputs => {
		GrothendieckWittClass => "beta" => {"the isomorphism class of a symmetric bilinear form"}
		},
	Outputs => { Ring => { "the base field of the Grothendieck-Witt class ", TT "beta" }},
	PARA {"Given the isomorphism class of a symmetric bilinear form, ", TT "beta", 
                ", this command outputs the base field of the form."},
	EXAMPLE lines ///
		 beta = gwClass(matrix(QQ,{{0,2},{2,0}}));
		 baseField beta
	 	 ///,
        }



document {
    Key => {(localAlgebraBasis, List, Ideal), localAlgebraBasis},
	Headline => "produces a basis for a local finitely generated algebra over a field k",
	Usage => "localAlgebraBasis(L,p)",
	Inputs => {
		List => "L" => {"list of polynomials ", TEX///$f=(f_1, \dots ,f_n)$///, " over the same ring"},
		Ideal => "p" => {"prime ideal of an isolated zero"}
	},
	Outputs => { List => {"a list of basis elements of the local k-algebra ", TEX///$Q_p(f)$/// }},
	PARA {"Given an endomorphism of affine space, ", TEX///$f=(f_1,\dots ,f_n)$///,
			", given as a list of polynomials called ", TT "L", " and the prime ideal of an isolated zero, this command returns a list of basis elements of the local k-algebra ", TEX///$Q_p(f)$///, "."},
	EXAMPLE lines ///
		 QQ[x,y];
		 f = {x^2+1-y,y};
		 p = ideal(x^2+1,y);
		 localAlgebraBasis(f,p) 
	 	 ///,
        }

document {
    Key => {(simplifyForm, GrothendieckWittClass), simplifyForm},
    Headline => "produces a simplified diagonal representative of a Grothendieck Witt class",
    Usage => "simplifyForm(beta)",
    Inputs => {
        GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over a field ", TEX///$k$///, "."},
    },
    Outputs => { 
	GrothendieckWittClass => {"a diagonal representative of the Grothendieck Witt class of the input form"},
	--String => {"The decomposition as a sum of hyperbolic and rank one forms."},
	},
    PARA {"Given a symmetric bilinear form ", TT"beta", " over a field ", TEX///$k$///, ", we diagonalize the form ", TT"beta", " and write it as a sum of some number of hyperbolic forms and some number of rank one forms."},
    EXAMPLE lines ///
    M = matrix(RR,{{2.091,2.728,6.747},{2.728,7.329,6.257},{6.747,6.257,0.294}});
    beta = gwClass(M);
    simplifyForm(beta)
    ///,
    PARA {"Over ", TEX///$\mathbb{R}$///, " there are only two square classes and a form is determined uniquely by its rank and signature [L05, II Proposition 3.2]. A form defined by the ", TEX///$3\times 3$///, " Gram matrix ", TT"M", " above is isomorphic to the form ", TEX///$\langle 1,-1,1\rangle $///, ", which decomposes as a sum ", TEX///$\mathbb{H}+\langle 1\rangle $///, "."},
    EXAMPLE lines ///
    M = matrix(GF(13),{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
    beta = gwClass(M);
    simplifyForm(beta)
    ///,
    PARA {"Over ", TEX///$\mathbb{F}_{q}$///, " forms can similarly be diagonalized and decomposed, in this case as ", TEX///$\mathbb{H} + \langle 1 \rangle + \langle -6 \rangle$///, "."},
    PARA{EM "Citations:"},
    UL{
	
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	},
}

document {
    Key => {(simplifyFormString, GrothendieckWittClass), simplifyFormString},
    Headline => "produces a simplified diagonal representative of a Grothendieck Witt class",
    Usage => "simplifyForm(beta)",
    Inputs => {
        GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over a field ", TEX///$k$///, "."},
    },
    Outputs => { 
	--GrothendieckWittClass => {"a diagonal representative of the Grothendieck Witt class of the input form"},
	String => {"The decomposition as a sum of hyperbolic and rank one forms."},
	},
    PARA {"Given a symmetric bilinear form ", TT"beta", " over a field ", TEX///$k$///, ", we diagonalize the form ", TT"beta", " and write it as a sum of some number of hyperbolic forms and some number of rank one forms."},
    EXAMPLE lines ///
    M = matrix(RR,{{2.091,2.728,6.747},{2.728,7.329,6.257},{6.747,6.257,0.294}});
    beta = gwClass(M);
    simplifyFormString(beta)
    ///,
    PARA {"Over ", TEX///$\mathbb{R}$///, " there are only two square classes and a form is determined uniquely by its rank and signature [L05, II Proposition 3.2]. A form defined by the ", TEX///$3\times 3$///, " Gram matrix ", TT"M", " above is isomorphic to the form ", TEX///$\langle 1,-1,1\rangle $///, ", which decomposes as a sum ", TEX///$\mathbb{H}+\langle 1\rangle $///, "."},
    EXAMPLE lines ///
    M = matrix(GF(13),{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
    beta = gwClass(M);
    simplifyFormString(beta)
    ///,
    PARA {"Over ", TEX///$\mathbb{F}_{q}$///, " forms can similarly be diagonalized and decomposed, in this case as ", TEX///$\mathbb{H} + \langle 1 \rangle + \langle -6 \rangle$///, "."},
    PARA{EM "Citations:"},
    UL{
	
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	},
}



document {
    Key => {(globalA1Degree, List), globalA1Degree},
    Headline => "Computes a global A1-Brouwer degree of a list of n polynomials in n variables over a field k",
    Usage => "globalA1Degree(L)",
    Inputs => {
	List => "L" => {"list of polynomials ", TEX///$f = (f_1, \ldots, f_n)$///, " in the polynomial ring ", TEX///$k[x_1,\ldots,x_n]$///, " over a field ", TEX///$k$///, "."}
	},
    Outputs => {
	GrothendieckWittClass => {"The class ", TEX///$\text{deg}^{\mathbb{A}^1}(f)$///, " in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."}
	},
    PARA{"Given an endomorphism of affine space ", TEX///$f=(f_1,\dots ,f_n) \colon \mathbb{A}^n_k \to \mathbb{A}^n_k$///, " with isolated zeros, we may compute its ", TEX///$\mathbb{A}^1$///, EM "-Brouwer degree", " valued in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."
	},
    PARA{"The ",
	TEX///$\mathbb{A}^1$///,
	EM "-Brouwer degree",
	"first defined by Morel [M12] is an algebrao-geometric enrichment of the classical topological Brouwer degree. Using the tools of motivic homotopy theory, one may associate to an endomorphism of affine space the isomorphism class of a symmetric bilinear form whose invariants encode geometric data about how the morphism transforms space."
        },
    PARA{"Such an association appears in the work of Eisenbud-Levine [EL77] and Khimshiashvili [K77], wherein the authors develop a symmetric bilinear form whose signature computes the local degree of a smooth map of real manifolds in the case where the Jacobian may vanish on an affine chart. This was proven to agree with Morel's ", TEX///$\mathbb{A}^1$///, "-Brouwer degree in work of Kass and Wickelgren [KW19]. A similar production of a symmetric bilinear form is given by work of Scheja and Storch [SS76], which develops a symmetric bilinear form attached to a complete intersection. This was also shown to align with the ", TEX///$\mathbb{A}^1$///, "-Brouwer degree in [BW23]."},
    PARA{"Following recent work of B. McKean and Pauli [BMP23], the ", TEX///$\mathbb{A}^1$///, "-Brouwer degree can be computed as a multivariate ", EM "Bezoutian bilinear form.", " The algorithms for producing such a form are developed here."},
    EXAMPLE lines ///
    QQ[x];
    f = {x^2+1};
    globalA1Degree(f)
    ///,
    PARA{"The previous example produces a rank two form with signature zero. This corresponds to the fact that the degree of the complex map ", TEX///$\mathbb{C}\to\mathbb{C},\ z\mapsto z^2$///, " has degree two, while the associated real map ", TEX///$\mathbb{R}\to\mathbb{R},\ x\mapsto x^2$///, " has global degree zero."},
    PARA{"Following [M21] we may think about the ",TEX///$\mathbb{A}^1$///, "-Brouwer degree ", TEX///$\text{deg}^{\mathbb{A}^1}(f)$///, " as a quadratically enriched intersection multiplicity of the hyperplanes ", TEX///$V(f_1)\cap \cdots \cap V(f_n).$///, " As a toy example, consider the curve ", TEX///$y=x(x-1)(x+1)$///, " intersecting the ", TEX///$x$///, "-axis."},
    EXAMPLE lines ///
    QQ[x,y];
    f = {x^3 - x^2 - y, y};
    globalA1Degree(f)
    ///,
    PARA{"The rank of this form is three, as cubics over the complex numbers have three roots counted with multiplicity. This form has signature one, which indicates that when the cubic intersects the ", TEX///$x$///, "-axis, when the three points of intersection are counted with a sign corresponding to a right hand rule, the sum equals one."},
    
    PARA{EM "Citations:"},
    UL{
	
	{"[BW23] T. Bachmann, K. Wickelgren, ", EM "Euler classes: six-functors formalism, dualities, integrality and linear subspaces of complete intersections,", " J. Inst. Math. Jussieu, 2023."},
	{"[EL77] D. Eisenbud, H. Levine, ", EM "An algebraic formula for the degree of a C^infinity map germ,", " Annals of Mathematics, 1977."},
	{"[K77] G. Khimshiashvili, ", EM "The local degree of a smooth mapping,", " Sakharth. SSR Mcn. Akad. Moambe, 1977."},
	{"[KW19] J. Kass, K. Wickelgren, ", EM "The class of Eisenbud-Khimshashvili-Levine is the local A1-Brouwer degree,", " Duke Math J., 2019."},
	{"[M21] S. McKean, ", EM "An arithmetic enrichment of Bezout's Theorem,", " Math. Ann., 2021."},
	{"[M12] F. Morel, ", EM "A1-Algebraic topology over a field,", " Springer Lecture Notes in Mathematics, 2012."},
	{"[BMP23] T. Brazelton, S. McKean, S. Pauli, ", EM "Bezoutians and the A1-Degree,", " Algebra & Number Theory, 2023."},
	{"[SS76] S. Scheja, S. Storch, ", EM "Uber Spurfunktionen bei vollstandigen Durchschnitten,", " J. Reine Angew. Math., 1975."},
	},
    
    }


document {
    Key => {(localA1Degree, List, Ideal), localA1Degree},
    Headline => "Computes a local A1-Brouwer degree of a list of n polynomials in n variables over a field k at a prime ideal in the zero locus",
    Usage => "locallA1Degree(L,p)",
    Inputs => {
	List => "L" => {"list of polynomials ", TEX///$f = (f_1, \ldots, f_n)$///, " in the polynomial ring ", TEX///$k[x_1,\ldots,x_n]$///, " over a field ", TEX///$k$///, "."},
	Ideal => "p" => {"a prime ideal ", TEX///$p \trianglelefteq k[x_1,\ldots,x_n]$///, " in the zero locus ", TEX///$V(f)$///, "."},
	},
    Outputs => {
	GrothendieckWittClass => {"The class ", TEX///$\text{deg}_p^{\mathbb{A}^1}(f)$///, " in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."}
	},
    PARA{"Given an endomorphism of affine space ", TEX///$f=(f_1,\dots ,f_n) \colon \mathbb{A}^n_k \to \mathbb{A}^n_k$///, " and an isolated zero ", TEX///$p\in V(f)$/// ,", we may compute its local ", TEX///$\mathbb{A}^1$///, EM "-Brouwer degree", " valued in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."
	},
    PARA{"For historical and mathematical background, see ", TO2(globalA1Degree, "global A1-degrees"), "."},
    EXAMPLE lines ///
    T1 = QQ[z_1..z_2];
    f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
    f1GD = globalA1Degree(f1);
    q=ideal {z_1,z_2};
    r=ideal {z_1-1,z_2^2-(9/85)};
    f1LDq= localA1Degree(f1,q)
    f1LDr= localA1Degree(f1,r)
    f1LDsum = gwAdd(f1LDq, f1LDr)
    ///,
    PARA{"The sum of the local A1-degrees is equal to the global A1-degree:"},
    PARA{"TODO --- add example after importing IsIsomorphic2"},
    
    }


document{
    Key => {(signature, GrothendieckWittClass), signature},
    Headline => "Outputs the signature of a symmetric bilinear form over the real or rational numbers",
    Usage => "signature(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"A symmetric bilinear form defined over ", TEX///$\mathbb{Q}$///, " or ", TEX///$\mathbb{R}$///, "."},
	},
    Outputs => {
	ZZ => "n" => {"The ", EM "signature", " of the symmetric bilinear form ", TEX///$\beta$///, "."},
	},
    PARA{"Given a symmetric bilinear form, after diagonalizing it, we can consider the number of positive entries minus the number of negative entries appearing along the diagonal. This is the ", EM "signature", " of a symmetric bilinear form, and is one of the primary invariants we use to classify forms. For more information see ", TO2(isIsomorphic2,"isIsomorphic2"), "."},
    EXAMPLE lines ///
    M = matrix(RR,{{0,0,1},{0,1,0},{1,0,0}});
    beta = gwClass(M);
    signature(beta)
    ///
    }


document{
    Key => {(isIsomorphic2, GrothendieckWittClass, GrothendieckWittClass), isIsomorphic2},
    Headline => "Determines whether two Grothendieck Witt classes are isomorphic over CC, RR, QQ, or a finite field.",
    Usage => "isIsomorphic2(alpha,beta)",
    Inputs => {
	GrothendieckWittClass => "alpha" => {"todo alpha"},
	GrothendieckWittClass => "beta" => {"todo beta"},
	},
    Outputs => {
	Boolean => {"returns true or false depending on whether two Grothendieck Witt classes are equal in the Grothendieck-Witt ring"},
	},
    PARA{"Given two matrices representing symmetric bilinear forms over a field ", TEX///$k$///, ", it is a fundamental question to ask when they are representing the same symmetric bilinear form, i.e. when they are equal in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///,"."},
    
    PARA{EM "Sylvester's Law of Inertia", " proves that any symmetric bilinear form can be diagonalized into a block sum of rank one symmetric bilinear forms. Since the rank one forms ", TEX///$\langle a \rangle \colon k \times k \to k$///, ", ", TEX///$(x,y) \mapsto axy$///, " and ", TEX///$\langle ab^2 \rangle \colon k \times k \to k$///, ", ", TEX///$(x,y) \mapsto ab^2xy$///, " differ by a change of basis in the ground field, it follows they are isomorphic (provided that ", TEX///$a,b\ne 0$///, "). Thus after diagonalizing a form, it suffices to consider the square class of each entry appearing along the diagonal. Consider the following example."},
    EXAMPLE lines ///
    alpha = gwClass(matrix(CC,{{2,3,1},{3,-1,0},{1,0,0}}))
    beta = gwClass(matrix(CC,{{2,4,-1},{4,5,7},{-1,7,9}}))
    isIsomorphic2(alpha,beta)
    ///,
    PARA{"The two forms are isomorphic since they can be diagonalized, after which they can be rewritten as the identity matrix after a chance of basis, since every nonzero element is a square class over ", TEX///$\mathbb{C}$///, " (the same is true for any quadratically closed field). Thus we have that the ", EM "rank", " of a form completely determines it over the complex numbers. That is, it provides an isomorphism ", TEX///$\text{GW}(\mathbb{C}) \to \mathbb{Z}$/// ,"."},
    PARA{"Over the reals, the story is a bit different. Since there are two classes over the reals, ", TEX///$ \mathbb{R}^\times / \left(\mathbb{R}^\times\right)^2 \cong \left\{\pm 1\right\}$///," we have a further invariant which classifies symmetric bilinear forms, called the ", TO2(signature, "signature"), ". This is computed as first diagonalizing, then taking the number of positive entries appearing on the diagonal minus the number of negative entries appearing on the diagonal."},
    EXAMPLE lines ///
    gamma = gwClass(matrix(RR,{{1,0,0},{0,-1,0},{0,0,1}}));
    signature(gamma)
    ///,
    PARA{"Rank and signature completely classify symmetric bilinear forms over the reals."},
    EXAMPLE lines ///
    delta = gwClass(matrix(RR,{{0,0,1},{0,1,0},{1,0,0}}));
    isIsomorphic2(gamma,delta)
    ///,
    PARA{"Over finite fields, rank is still an invariant of a form, however signature no longer makes sense as the field is not totally ordered. Instead we consider the ", EM "discriminant", " of the non-degenerate symmetric bilinear form, which is the determinant of any Gram matrix representing the form. The discriminant is well-defined once we consider its target as landing in square classes of the field. ", TEX///$\text{GW}(\mathbb{F}_q) \to \mathbb{F}_q^\times / \left(\mathbb{F}_q^\times\right)^2$///, ". Over finite fields we look to compare both rank and discriminant of a form."},
    EXAMPLE lines///
    alphaF = gwClass(matrix(GF(7),{{1,2,2},{2,0,1},{2,1,5}}))
    betaF = gwClass(matrix(GF(7),{{2,5,1},{5,6,1},{1,1,3}}))
    gammaF = gwClass(matrix(GF(7),{{0,2,4},{2,3,3},{4,3,1}}))
    det(alphaF.matrix)    
    det(betaF.matrix)
    det(gammaF.matrix)
    isIsomorphic2(alphaF,betaF)
    isIsomorphic2(alphaF,gammaF)
    isIsomorphic2(betaF,gammaF)
    legendreBoolean(det(betaF.matrix)) ==legendreBoolean(det(gammaF.matrix))
    ///,
    PARA{"ok figure out the error above"},
    PARA{"Over the rationals, further invariants must be considered. We first check if the rank, discriminant, and signature (when considered as a real form) all agree. If so, we must further check whether the ", EM "Hasse-Witt invariants", " agree at all primes. This is an instance of the ", EM "Hasse-Minkowski principle", " which states that quadratic forms are isomorphic over a global field if and they are isomorphic over all its completions (see [S73, III Theorem 7] or [L05, VI.3.3])."},
    PARA{"The ", EM "Hasse-Witt invariant", " of a diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " over a field ", TEX///$K$///, " is defined to be the product ", TEX///$\prod_{i<j} \left( \phi(a_i,a_j) \right)$///, " where ", TEX///$\phi \colon K \times K \to \left\{\pm 1\right\}$///, " is any ", EM "symbol", " (see e.g. [MH73, III.5.4] for a definition). It is a classical result of Hilbert that over a local field of characteristic not equal to two, there is one and only symbol, ", TEX///$(-,-)_p$///,  " called the ", EM "Hilbert symbol", " ([S73, Chapter III]) computed as follows:"},
    PARA{TEX///$(a,b)_p = \begin{cases} 1 & z^2 = ax^2 + by^2 \text{ has a nonzero solution in } K^3 \\ -1 & \text{otherwise.} \end{cases}$///},
    PARA{"Consider the following example, where we observe that ", TEX///$z^2 = 2x^2 + y^2$///," does admit nonzero solutions mod 7, in particular ", TEX///$(x,y,z) = (1,0,3)$///, ":"},
    EXAMPLE lines///
    hilbertSymbol(2,1,7)
    ///,
    PARA{"The Hasse invariant will be 1 for almost all primes. In particular after diagonalizing a form ", TEX///$\beta \cong \left\langle a_1,\ldots,a_n\right\rangle$///, " then the Hasse invariant at a prime ", TEX///$p$///, " will be 1 automatically if ", TEX///$p\nmid a_i$///, " for all ", TEX///$i$///, ". Thus we only have finitely many Hasse invariants to compare for any pair of symmetric bilinear forms."},
    EXAMPLE lines///
    alphaQ = gwClass(matrix(QQ,{{1,4,7},{4,3,2},{7,2,-1}}))
    betaQ = gwClass(matrix(QQ,{{0,0,1},{0,2,7},{1,7,3}}))
    isIsomorphic2(alphaQ,betaQ)
    ///,
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	{"[MH73] Milnor and Husemoller, ", EM "Symmetric bilinear forms,", " Springer-Verlag, 1973."},
    },
}




document{
    Key => {(hilbertSymbol, ZZ,ZZ,ZZ), hilbertSymbol},
    Headline => "Computes the Hilbert symbol of two integers at a prime",
    Usage => "hilbertSymbol(a,b,p)",
    Inputs => {
	ZZ => "a" => {"Any integer, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	ZZ => "b" => {"Any integer, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	ZZ => "p" => {"Any prime number."},
	},
    Outputs => {
	ZZ => "(a,b)_p" => {"The ", EM "Hilbert symbol ", TEX///$(a,b)_p$///, "."},
	},
    PARA{"The ", EM "Hasse-Witt invariant", " of a diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " over a field ", TEX///$K$///, " is defined to be the product ", TEX///$\prod_{i<j}  \phi(a_i,a_j)$///, " where ", TEX///$\phi \colon K \times K \to \left\{\pm 1\right\}$///, " is any ", EM "symbol", " (see e.g. [MH73, III.5.4] for a definition). It is a classical result of Hilbert that over a local field of characteristic not equal to two, there is one and only symbol, ", TEX///$(-,-)_p$///,  " called the ", EM "Hilbert symbol", " ([S73, Chapter III]) computed as follows:"},
    PARA{TEX///$(a,b)_p = \begin{cases} 1 & z^2 = ax^2 + by^2 \text{ has a nonzero solution in } K^3 \\ -1 & \text{otherwise.} \end{cases}$///},
    PARA{"Consider the following example, where we observe that ", TEX///$z^2 = 2x^2 + y^2$///," does admit nonzero solutions mod 7, in particular ", TEX///$(x,y,z) = (1,0,3)$///, ":"},
    EXAMPLE lines///
    hilbertSymbol(2,1,7)
    ///,
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[MH73] Milnor and Husemoller, ", EM "Symmetric bilinear forms,", " Springer-Verlag, 1973."},
    },
}



document {
    Key => {(diagonalForm, GrothendieckWittClass), diagonalForm},
	Headline => "produces a diagonalized form for any Grothendieck-Witt class",
	Usage => "diagonalForm(beta)",
	Inputs => {
		GrothendieckWittClass => "beta" => {"any class in ", TEX///$\text{GW}(k)$///," where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."}
	},
	Outputs => { GrothendieckWittClass => {"a form isomorphic to ", TEX///$\beta$///, " with a diagonal Gram matrix"}},
	PARA {"test"},
	EXAMPLE lines ///
	beta = gwClass(matrix(QQ,{{0,0,2},{0,2,0},{2,0,0}}));
	diagonalForm(beta)
        ///,
	PARA{"Note that the  ", TO2(GrothendieckWittClass, "GrothendieckWittClass"), " type caches diagonal versions of a form once they've been computed. We can recover this quickly in the following way."},
	EXAMPLE lines///
	beta.cache.diagonalForm
	///,
	SeeAlso => {"diagonalize"}
	}




document{
    Key => {(isIsotropic, GrothendieckWittClass), isIsotropic},
    Headline => "Determines whether a Grothendieck-Witt class is isotropic",
    Usage => "isIsotropic(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(k)$///, " where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."},
	},
    Outputs => {
        Boolean => {"Whether ", TEX///$\beta$///, " is isotropic"},
	},
    PARA{"Recall a symmetric bilinear form ", TEX///$\beta$///, " is said to be ", EM "isotropic", " if there exists a nonzero vector ", TEX///$v$///, " for which ", TEX///$\beta(v,v) = 0$///, ". Witt's decomposition theorem implies that a non-degenerate symmetric bilinear form decomposes uniquely into an isotropic and an anisotropic part. Certifying (an)isotropy is then an important computational problem when working with the Grothendieck-Witt ring."},
    PARA{"Over ", TEX///$\mathbb{C}$///, ", any form of rank two or higher contains a copy of the hyperbolic form, and hence is isotropic. Thus we can determine isotropy simply by a consideration of rank."},
    EXAMPLE lines///
    isIsotropic(gwClass(matrix(CC,{{3}})))
    isIsotropic(gwClass(matrix(CC,{{2,0},{0,5}})))
    ///,
    PARA{"Forms over ", TEX///$\mathbb{R}$///, " are anisotropic if and only if all its diagonal entries are positive or are negative."},
    EXAMPLE lines///
    isIsotropic(gwClass(matrix(RR,{{3,0,0},{0,5,0},{0,0,7}})))
    isIsotropic(gwClass(matrix(RR,{{0,2},{2,0}})))
    ///,
    PARA{"Over finite fields, a form is anisotropic so long as it is nondegenerate, of rank ", TEX///$\le 2$///," and not isomorphic to the hyperbolic form."},
    EXAMPLE lines///
    isIsotropic(gwClass(matrix(GF(7),{{1,0,0},{0,1,0},{0,0,1}})))
    isIsotropic(gwClass(matrix(GF(7),{{3,0},{0,3}})))
    ///,
    PARA{"Over ", TEX///$\mathbb{Q}$///, " things become a bit more complicated. We can exploit the local-to-global principle for isotropy (the ", EM "Hasse-Minkowski principle", "), which states that a form is isotropic over ", TEX///$\mathbb{Q}$///, " if and only if it is isotropic over all its completions, meaning all the ", TEX///$p$///, "-adic numbers and ", TEX///$\mathbb{R}$///, " [L05, VI.3.1]. We note, however, the classical result that all forms of rank ", TEX///$\ge 5$///, " in ", TEX///$\mathbb{Q}_p$///, " are isotropic [S73, III Theorem 6]. Thus isotropy in this range of ranks is equivalent to checking it over the real numbers."},
    EXAMPLE lines///
    beta = gwClass(matrix(QQ,{{1, 0, 2, 0, 3}, {0, 6, 1, 1, -1},{2, 1, 5, 2, 0}, {0, 1, 2, 4, -1}, {3, -1, 0,-1, 1}}));
    isIsotropic(beta)
    diagonalForm(beta)
    ///,
    PARA{"For forms of rank ", TEX///$\le 4$///, " there are simple criteria for isotropy over ", TEX///$\mathbb{Q}_p$///, " which can be found in [S73, III Theorem 6]. As an example for rank 3 forms, isotropy of a form ", TEX///$\beta \in \text{GW}(\mathbb{Q})$///, " over ", TEX///$\mathbb{Q}_p$///," is equivalent to the statement that ", TEX///$(-1,-\text{disc}(\beta))_p = H(\beta)$///, " where ", TEX///$H(\beta)$///, " denotes the Hasse-Witt invariant attached to ", TEX///$\beta$///, " and ", TEX///$(-,-)_p$///," is the ", TO2(hilbertSymbol, "Hilbert Symbol"), "."},
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
    },
    
}


document{
    Key => {(isAnisotropic, GrothendieckWittClass), isAnisotropic},
    Headline => "Determines whether a Grothendieck-Witt class is anisotropic",
    Usage => "isAnisotropic(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(k)$///, " where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."},
	},
    Outputs => {
        Boolean => {"Whether ", TEX///$\beta$///, " is anisotropic"},
	},
    PARA{"This is the negation of the boolean-valued ", TO2(isIsotropic,"isIsotropic"), ". See documentation there."},
    
}

------------------
-- TESTING
------------------

-- Diagonal form testing
TEST ///
print("diagonal form testing");
M1=matrix(RR, {{0, 1}, {1, 0}});
G1=gwClass(M1);
M2=diagonalForm(G1);
assert(M2.matrix===matrix(RR, {{1, 0}, {0, -1}}));
///

TEST ///
M3=matrix(CC, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
assert(M4.matrix===matrix(CC, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
///

TEST ///
M5=matrix(RR, {{1, 2, 3, 4}, {2, 4, 5, 16}, {3, 5, 7, 8}, {4, 16, 8, 19}});
G3=gwClass(M5);
M6=diagonalForm(G3);
A=matrix(RR, {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}});
assert(M6.matrix===A);
///

TEST ///
M3=matrix(QQ, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
assert(M4.matrix===matrix(QQ,{{1, 0, 0}, {0, -4, 0}, {0, 0, 1/4}}));
///

-- gwTypeTest.m2
TEST ///
M = matrix(QQ,{{1,0},{0,1}});
N = matrix(QQ, {{1, 2}, {3, 4}})
beta = gwClass(M);
gamma = gwClass(N);
assert(baseField(beta) === QQ)
assert(beta.matrix === M)
--Operations within GW-classes
A = gwAdd(beta, gamma);
B = gwMultiply(beta, gamma);
assert(A.matrix === matrix(QQ, {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 2}, {0, 0, 3, 4}}));
assert(B.matrix === matrix(QQ, {{1, 2, 0, 0}, {3, 4, 0, 0}, {0, 0, 1, 2}, {0, 0, 3, 4}}));
---non well-defined GW-classes
M'=matrix(ZZ, {{1, 0}, {0, 1}});
N'=matrix(QQ, {{1, 1}, {1, 1}});
theta=gwClass(M');
sigma=gwClass(N');
assert(isWellDefined(theta) === false);
assert(isWellDefined(sigma) === false);
///


-- degreeTesting.m2
TEST ///
T1 = QQ[x]
f = {x^2}
beta = globalA1Degree(f)

gamma = gwClass(matrix(QQ,{0,1},{1,0}))
assert(beta==gamma)
///


TEST ///
T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
f1GD = globalA1Degree(f1);
f1GDmat = f1GD.matrix;
assert((wittDecomp(f1GDmat))==(3,matrix(QQ,{{}})));
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
f1LDq= localA1Degree(f1,q);
f1LDr= localA1Degree(f1,r);
f1LDsum = gwAdd(f1LDq, f1LDr);
assert(isIsomorphic2(f1LDsum, f1GD));



T2 = QQ[w];
f2 = {w^4 + w^3 - w^2 - w};
f2GD= globalA1Degree(f2);
f2GDmat = f2GD.matrix;
assert(wittDecomp(f2GDmat)==(2,matrix(QQ,{{}})));

p=ideal {w+1};
f2LDp = localA1Degree(f2, p);
f2LDpmat = f2LDp.matrix;
assert(wittDecomp(f2LDpmat)==(1,matrix(QQ,{{}})));
s=ideal{w-1};
f2LDs = localA1Degree(f2, s);
t=ideal{w};
f2LDt = localA1Degree(f2, t);
f2LDsum = gwAdd(gwAdd(f2LDp, f2LDs),f2LDt);
assert(isIsomorphic2(f2LDsum, f2GD));
///
    


end
