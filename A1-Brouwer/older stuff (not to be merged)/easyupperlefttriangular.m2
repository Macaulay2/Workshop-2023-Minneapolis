--Joel Louwsma

--return true if can easily detect that square matrix is upper left triangular after permuting basis vectors
--otherwise return false.

easyUpperLeftTriangular = method()

easyUpperLeftTriangular (Matrix) := (Boolean) => (A) -> (
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    while numRows(A) >= 2 do (
	goodRow = numRows(A);
        zeroRow = numRows(A);
        for i from 0 to (numRows(A)-1) do (
            zeroCount = 0;
    	    for j from 0 to (numRows(A)-1) do (
	        if A_(i,j) == 0 then (
	            zeroCount = zeroCount+1;
		    );
		);
	    if zeroCount == (numRows(A)) then (
		zeroRow = i;
		break;
		);            
	    if (zeroCount == (numRows(A)-1) and A_(i,i) == 0) then (
		goodRow = i;
		break;
		);
	    );
	--if no row has numRows(A)-1 zeros, then need more complicated analysis than this function does
        if (goodRow != numRows(A)) then (
	    for j from 0 to (numRows(A)-1) do (
	        if A_(goodRow,j) != 0 then (
		    goodColumn = j;
		    break;
		    );
	        );
	    --if we have a row with numRows(A)-1 zeros, remove this row and column and the corresponding row and column and continue
	    A = submatrix'(A,{goodRow,goodColumn},{goodRow,goodColumn});
            )
        --if we have a row with numRows(A) zeros, remove this row and column and continue
        else if (zeroRow != numRows(A)) then (
            A = submatrix'(A,{zeroRow},{zeroRow});
            )
	else (
	    return false;
	    );
        );
     --if we get this far, numRows(A) is 0 or 1, meaning after a permutation of basis vectors the orginal matrix is upper left triangular
     true
     )
