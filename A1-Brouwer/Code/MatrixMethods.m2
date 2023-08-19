-----------------------
-- Matrix manipulations
-----------------------

-- Input: Two matrices.
-- Output: The block sum of possibly empty matrices.

-- Note: The block sum "++" of a zero matrix with something else outputs the wrong thing.
-- This internal method allows us to call block sums of possibly empty matrices.

safeBlockSum = method()
safeBlockSum (Matrix, Matrix) := Matrix => (A,B) -> (
    if numColumns A == 0 then return B;
    if numColumns B == 0 then return A;    
    return A ++ B 
    )

-- Input: A matrix.
-- Output: True if the matrix is a square; false otherwise.

isSquare = method()
isSquare (Matrix) := Boolean => M -> (
    numRows(M) == numColumns(M)
    )

-- Input: A matrix
-- Output: True if the matrix is a square and symmetric; false otherwise.

isSquareAndSymmetric = method()
isSquareAndSymmetric (Matrix) := Boolean => M -> (
    transpose(M) == M
    )

-- Input: A matrix.
-- Output: True if the matrix represents a degenerate bilinear form; false otherwise.

isDegenerate = method()
isDegenerate (Matrix) := Boolean => M ->(
    det(M) == 0
    )

-- Input: A matrix.
-- Output: True if the matrix is upper-left triangular, meaning that all the entries below the main antidiagonal are zero; false otherwise.

isUpperLeftTriangular = method()
isUpperLeftTriangular (Matrix) := Boolean => M -> (
    if not isSquare(M) then error "Error: matrix isn't square";
    n := numRows(M);
    for i from 1 to n - 1 do(
	for j from 0 to i - 1 do(
    	-- If any entry in this range is nonzero then the matrix isn't upper left triangular
		if not M_(i,j) == 0 then(
		    return false
		    );
		
        );
    );
    true
    )
 
-- Input: A square matrix.
-- Output: True if the matrix is diagonal; false otherwise. 

isDiagonal = method()
isDiagonal (Matrix) := Boolean => M -> (

    if not isSquare(M) then error "Error: matrix is not a square";

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

-- Input: A symmetric matrix.
-- Output: A diagonal matrix congruent to the original matrix.

congruenceDiagonalize = method()
congruenceDiagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    k := ring AnonMut;
    if isField k == false then error "Error: expected matrix entries from a field";
    if not isSquareAndSymmetric(AnonMut) then error "matrix is not symmetric";
    
    -- If the matrix is already diagonal then return it
    if isDiagonal(AnonMut) == true then(
	return AnonMut
	);
    
    -- Otherwise we iterate through columns and rows under the diagonal, and perform row operations followed by the corresponding
    -- transpose operation on columns in order to reduce to a diagonal matrix congruent to the original
    A := mutableMatrix AnonMut;
    n := numRows(A);
    for col from 0 to (n - 1) do (
	-- If diagonal entry in column "col" is zero
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
                -- Better alternative for keeping entries smaller in A in case A_(row,row)!=0
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
                temp:=A_(row,col);
                -- More row reduction make every entry below A_(col,col) is zero
                rowAdd(A,row,-temp/A_(col,col),col);
	        -- Column reduction to keep reduced matrix congruent
                columnAdd(A,row,-temp/A_(col,col),col);
                );
            );
        );
    return matrix A 
    )
