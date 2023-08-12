

-- Input: A square matrix.
-- Output: The radical of a quadratic space.

-- Note: This is reliant on the congruenceDiagonalize() method being applicable for singular matrices.

truncateRadical = method()
truncateRadical(Matrix) := (Matrix) => (A) -> (
    truncatedMatrix := mutableMatrix A;
    if not isSquare(A) then error ("Input is not a square matrix");
   
    truncatedMatrix = mutableMatrix congruenceDiagonalize(A);
    foundRadical := false;
    for i from 0 to (numRows(A) - 1) do (
        if truncatedMatrix_(i, i) == 0 then (
            foundRadical = true;
            break
            );
        error ("The quadratic space does not have a radical!");
        );
    if foundRadical === true then (
        n:=numRows(A) - 1;
        for i from 0 to n do (
            truncatedMatrix = mutableMatrix submatrix'(matrix truncatedMatrix, {i}, {i});
            if (n > 0) then (n = n - 1;)
            else (break);
            );
        B := matrix truncatedMatrix;
        return B;
        );
    )

-- Input: A diagonal matrix.
-- Output: The number of times we split off any hyperbolic forms < a > + < - a > as well as the smaller matrix with none of them.

splitOffObviousHyperbolic = method()
splitOffObviousHyperbolic (Matrix) := (ZZ,Matrix) => (A) -> (
    
    -- Matrix must be symmetric
    if not isSquareAndSymmetric(A) then error "Matrix is not symmetric";

    foundHyperbolic := 0;
    remainingMatrix := A;
    for i from 0 to (numRows(A) - 1) do (
        for j from (i + 1) to (numRows(A) - 1) do (
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
    numberHyperbolics := 0;
    notFinished := 1;
    while (notFinished == 1) do (
        currentState := splitOffObviousHyperbolic(A);
        notFinished = currentState_0;
        numberHyperbolics = numberHyperbolics + notFinished;
        A = currentState_1;
        );
    return (numberHyperbolics,A);
    )

-- Input: A symmetric matrix with rational entries.
-- Output: The number of times we split off any hyperbolic forms < a > + < - a > as well as the smaller matrix with none of them.

-- Note: Takes in symmetric matrix over QQ and diagonalizes, removes squares from entries, and splits off hyperbolic forms that immediately appear as < a > + < - a >
rationalSimplify = method()
rationalSimplify (Matrix) := (ZZ,Matrix) => (A) -> (
    -- Matrix must be symmetric
    if not isSquareAndSymmetric(A) then error "Matrix is not symmetric";

    -- congruenceDiagonalize the matrix
    B := mutableMatrix(congruenceDiagonalize(A));
    -- Replace entry with smallest magnitude integer in square class
    for i from 0 to (numRows(B)-1) do (
        B_(i,i) = squarefreePart(B_(i,i));
        );
    C := matrix B;
    -- Split off hyperbolic forms < a > + < - a >
    return(splitOffObviousHyperbolics(C))
    )



-- Nikita Borisov and Frenly Espino

-- Given a matrix A, WittDecomp decomposes A as a sum nH + Q, where n = number of hyperbolic forms, and Q is (presumed to be) anisotropic.  

-- Input:  A square symmetric matrix A over an exact field (not RR or CC)
-- Output: (n, Q), n = integer giving number of hyperbolic spaces, Q=anisotropic forms.

-- Future Work: Should be able to input a bound.

WittDecomp = method()
WittDecomp (Matrix) := (ZZ,Matrix) => (A) -> (
    k := ring A;   
  
    -- Add error in case the base field is RR or CC
    if (instance(k,InexactFieldFamily) or instance(k,RealField) or instance(k,ComplexField)) then error "Error: base field is inexact, use WittDecompInexact() instead";
    
    -- Rank of matrix
    n := numRows(A);
    x := symbol x;
    -- Local variable ring
    R := k[x_0..x_(n - 1)];
    
    -- Create quadratic from f from matrix A
    f := sum (
        for i from 0 to (n - 1) list (
            sum (for j from 0 to (n - 1) list (A_(i,j)*x_i*x_j))
            )
        );
 
    -- will use rationalPoints package to see if f has a zero.  If so, A is isotropic.  
    -- rationalPoints package seeks a zero upto variable bound
    
    use k;
    solnPt := new MutableList;
    solnFound := false;
    for bound from 1 to 10 do (
        solns := new MutableList from rationalPoints(ideal(f),Bound => bound);
        if (#solns > 1) then(
            if solns#0 == toList(n:0) then (solnPt = solns#1) else (solnPt = solns#0);
            solnFound = true;
            break;
        );
    );


    -- If no solutions found, then we assume the form is anisotropic
    if (not solnFound) then (return (0,A));
    
    -- If solution is found for a rank 2 form, then the form is purely hyperbolic
    if ((n == 2) and  (det(A) == 0))then (error "Matrix singular, run WittDecompGeneral instead" );
    if (n == 2) then (return (1,matrix(k,{{}})));

    -- If found a solution, record it as row matrix z. Then find vector y such that < z,y > <>0.  
    -- z and y then generate a hyberbolic plane. 
    -- z as a row matrix
    z := matrix{toList solnPt};
    zA := z*A;
    y := new MutableMatrix from matrix{toList(n:(0/1))};
    for i from 0 to (n - 1) do (
        if (zA_(0,i) != 0) then (y_(0,i) = 1; break;);
    );
    
    -- Now z and y span a copy of |H in the bilinear form
    -- We need to find a basis of vectors orthogonal (wrt bilinear form) to x and y 
    -- An (n - 2) x n matrix.
    orthoComp := gens kernel((z||matrix(y))*A);
    
    -- Now recursively apply WittDecomp to orthoComp^T*A*orthoComp a (n - 2) x (n - 2) Gram matrix
    subComputation := WittDecomp(transpose(orthoComp)*A*orthoComp);
    
    -- subComputation_0 gives number of hyperbolic forms in (n - 2) x (n - 2) subform  
    -- 1+ subComputation_0 is the number of hyperbolic forms in A
    -- subComputation_1 is the anisotropic part of A and (also) the subform part
    
    return (1+subComputation_0, subComputation_1);
    )


-- WittDecomp method for InexactFieldFamily

-- WittDecompInexact calculates the Witt decomposition for matrix A over the fields kk = RR, CC
-- Function assumes that A has maximal rank 

WittDecompInexact = method()
WittDecompInexact (Matrix) := (ZZ,Matrix) => (A) -> (
    k := ring A;
    
    -- Checks that we have RR or CC as our field
    if not (instance(k,RealField) or instance(k,ComplexField) or instance(k,InexactFieldFamily)) then error "Error: base field is not RR or CC";
    
    n := numRows(A); --rank of matrix
    
    -- If k is the complex numbers, witt decomposition depends only on rank
    -- If rank is even, then matrix decomposes into n/2 hyberbolic forms with no anisotropic parts
    -- If rank is odd, matrix decomposes into (n-1)/2 hyperbolic forms with 1 x 1 anisotropic part
    if (k === CC or instance(k,ComplexField)) then (
        if (n%2 == 0) then(return (n//2,matrix(CC,{{}})))
        else return (n//2,id_(k^1));
        );
    
    -- If k is the real numbers, witt decomposition depends on rank and signature
    if (k === RR or instance(k,RealField)) then (
        diagA := congruenceDiagonalize(A);
	-- for loop counts the number of positive diagonal entries of diagA
        posEntries := 0;
	-- for loop counts the number of negative diagonal entries
        negEntries := 0;
	for i from 0 to (n - 1) do(
            if diagA_(i,i) > 0 then(
                posEntries = posEntries+1;
            );
	    if diagA_(i,i) < 0 then(
                negEntries = negEntries+1;
            );
        );

        if (posEntries + negEntries > n) then (error "A is singular");
	-- Witt index is given by how many positive-negative diagonal entry pairs exist
        wittIndex := min(posEntries,negEntries);
        signature := posEntries - negEntries; 
        if signature == 0 then (return (wittIndex,matrix(RR,{{}})))
	-- signature characterizes anisotropic part
        else if signature > 0 then ( return (wittIndex, id_(k^(signature))))
        else return (wittIndex, -id_(k^(-signature)));
        );
    );

