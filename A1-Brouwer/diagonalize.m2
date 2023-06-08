--many people, essentially everyone

--diagonalize method
--given a symmetric matrix, this function outputs a diagonal matrix congruent to original matrix
diagonalize = method()

diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
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
