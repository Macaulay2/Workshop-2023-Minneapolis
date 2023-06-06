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

