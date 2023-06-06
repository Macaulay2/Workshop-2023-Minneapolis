--given a symmetric matrix, find a congruent diagonal matrix
diagonalize = method()

diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    if AnonMut != transpose(AnonMut) then (
        error "Matrix is not symmetric";
	);
    A := mutableMatrix AnonMut;
    n:=numRows(A);
    for col from 0 to (n-1) do (
        if A_(col,col) == 0 then (
            for row from col+1 to n-1 do (
                --try to find a nonzero entry in the column and use it to make the diagonal entry nonzero
                if A_(row,col) != 0 then (
                    rowAdd(A,col,1,row);
                    columnAdd(A,col,1,row);
                    break;
                );
            );
        );
        --make all entries in column other than diagonal entry 0
        if A_(col,col) != 0 then (
            for row from (col+1) to (n-1) do (
                temp:=A_(row,col);
                rowAdd(A,row,-temp/A_(col,col),col);
                columnAdd(A,row,-temp/A_(col,col),col);
            );
        );
    );
    return matrix A
)
