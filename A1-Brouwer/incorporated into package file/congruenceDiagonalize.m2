--Authors: many people, essentially everyone

--congruenceDiagonalize method
--Given a symmetric matrix, this function outputs a diagonal matrix congruent to original matrix
congruenceDiagonalize = method()

congruenceDiagonalize (Matrix) := (Matrix) => (AnonMut) -> (
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


