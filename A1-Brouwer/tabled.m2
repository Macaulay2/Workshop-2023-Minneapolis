
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
