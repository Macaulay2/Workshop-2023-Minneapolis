
-- Input: A symmetric matrix.
-- Output: A diagonal matrix congruent to the original matrix, capable of diagonalizing over rings (algorithm has no divisions).

congruenceDiagonalizeOverInt = method()
congruenceDiagonalizeOverInt (Matrix) := (Matrix) => (AnonMut) -> (
    if not isSquareAndSymmetric(AnonMut) then error "matrix is not symmetric";
    
    A := mutableMatrix AnonMut;
    n := numRows(A);
    for col from 0 to (n - 1) do (
        if A_(col,col) == 0 then (
            for row from col + 1 to n - 1 do ( 
		-- Scan for nonzero entries below the diagonal entry
                if A_(row,col) != 0 then (
                    if A_(row,row) == 0 then (
		        -- Row reduction to make A_(col,col) nonzero
                        rowAdd(A,col,1,row);
		        -- Column reduction to keep reduced matrix congruent to original matrix
                        columnAdd(A,col,1,row);
                        )
                    else (
		        -- Row and column swaps to make A_(col,col) nonzero
                        rowSwap(A,col,row);
                        columnSwap(A,col,row);
                        );
                    break;
                    );
                );
            );
	
        -- Now A_(col,col) != 0 unless there was a zero row/column and we use it to clear the column below
        if A_(col,col) != 0 then (
            for row from (col+1) to (n-1) do (
                temp := A_(row,col);
		
		-- Multiply row row by A_(col,col)
                rowMult(A,row,A_(col,col)); 
		
		-- Column multiplication to keep reduced matrix congruent
                columnMult(A,row,A_(col,col));
		
		-- More row reduction make every entry below A_(col,col) is zero
                rowAdd(A,row,-temp,col); 
		
		-- Column reduction to keep reduced matrix congruent
                columnAdd(A,row,-temp,col); 
                );
            );
        );
    return matrix A 
    )




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




-- TODO can we delete?
    
-- Input: (Q): Quadratic form Q given by list of diagonal elements.  
--      Assume list constists of integers
-- Output:  (Rank, Disc, Signature, Hasse Invariant for all primes p when not 1)    
invariantFormQ =method()
invariantFormQ (List):= (ZZ, ZZ, List, List) => (f) -> (
    -- currently will export the discriminant as a square free integer
    -- Note:  Still need a way to treat two integers as defining the same discriminant, if they differ by a
    --  square in Z_p.  They need to have the same parity of power of prime p, and the quotient of their 
    -- (prime-to-p) parts must define a unit square in Z_p^*.

    len:=#f;
    for i from 0 to (len-1) do (
    if (f_i==0) then (error "Error: Form is degenerate");
    if not liftable(f_i,ZZ ) then (error "Error: Diagonal elements of form should be integers");
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
    if ((not H#?2) and HasseWittInvariant(f, 2) == -1) then l=append(l, 2);
    for i from 0 to #k-1 do (
    if  (HasseWittInvariant(f, k_i) == -1) then (
        l=append(l,k_i);
        );
    );
    l=sort(l);
    return (a, b, c, l);
    );