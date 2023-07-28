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
    "globalA1Degree",
    "localA1Degree",
    "localAlgebraBasis"
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
        H:=hashTable(factor(numerator(n)*denominator(n)));
        return product(apply(keys(H),p->p^(H#p%2)))
        );
    if n < 0 then (
        H:=hashTable(factor(numerator(-n)*denominator(-n)));
        return -product(apply(keys(H),p->p^(H#p%2)))
        );
    )

squarefreePart (ZZ) := (ZZ) => (n) -> (
    if n==0 then (
        return 0
        );
    if n > 0 then (
        H:=hashTable(factor(n));
        return product(apply(keys(H),p->p^(H#p%2)))
        );
    if n < 0 then (
        H:=hashTable(factor(-n));
        return -product(apply(keys(H),p->p^(H#p%2)))
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
		    --Row reduction to make A_(col,col) nonzero
                    rowAdd(A,col,1,row);
		    --Column reduction to keep reduced matrix congruent to original matrix
                    columnAdd(A,col,1,row);
                    break;
                );
            );
        );
        --If nonzero entry at or below A_(col,col) is found, we use it to clear the column below
        if (A_(col,col)!=0) then (
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
        if A_(col,col) == 0 then ( --if diagonal entry in column "col" is zero
            for row from col+1 to n-1 do ( 
                if A_(row,col) != 0 then ( --scan for nonzero entries below the diagonal entry
                    rowAdd(A,col,1,row); --row reduction to make A_(col,col) non-zero
                    columnAdd(A,col,1,row); --column reduction to keep reduced matrix congruent to original matrix
                    break;
                );
            );
        );
        --if non-zero entry at or below A_(col,col) was found we use it to clear the column below
        if (A_(col,col)!=0) then (
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
        notFinished:= currentState_0;
        numberHyperbolics:= numberHyperbolics + notFinished;
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

-- Future Work: Shoudl be able to input a bound.

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

simplifyForm = method()
simplifyForm (GrothendieckWittClass) := (GrothendieckWittClass, String) => beta -> (
    print("got here");
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
            numHypFormsCCOdd := floor n/2;
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



-- New type of hash table called "GrothendieckWittType"
-- GW-type.m2


---------------
-- Operations with matrices
---------------

-- Diagonalization of symmetric bilinear forms
-- diagonalize.m2


-- basic booleans about matrices
-- matrixBooleans.m2


-- Checks if a Gram matrix can easily be seen to be upper left triangular
-- easyUpperTriangular.m2

---------------
-- Operations with k-algebras
---------------

-- localAlgebraBasis.m2


---------------
-- Diagonalization over QQ
---------------

-- Takes in a rational number or integer and outputs the smallest magnitude integer in its square class
-- squarefreepart.m2

-- Given diagonal matrix, split off any <a>+<-a> and return number of times we can do this as well as smaller matrix with none of these
-- splitoffobvioushyperbolics.m2

-- Takes in symmetric matrix over QQ and diagonalizes, removes squares from entries, and splits off hyperbolic forms that immediately appear as <a> + <-a>
-- rationalsimplify.m2





-- Checks if 
-- easyIsomorphicGW




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
	 	 ///
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
