needsPackage "RationalPoints2"

load "./GW-type.m2"

--diagonalize method
diagonalize = method()

diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    A := mutableMatrix AnonMut;
    n:=numRows(A);
    for col from 0 to (n-1) do (
        if A_(col,col) == 0 then (
            for row from col+1 to n-1 do (
                if A_(row,col) != 0 then (
                    --we have found non-zero entry
                    rowAdd(A,col,1,row);
                    columnAdd(A,col,1,row);
                    break;
                );
            );
            if A_(col,col)==0 then (print "Error: Matrix was singular"; return A;);
        );
        --entry in A_(col,col) is non-zero at this point
         for row from (col+1) to (n-1) do (
            temp:=A_(row,col);
            rowAdd(A,row,-temp/A_(col,col),col);
            columnAdd(A,row,-temp/A_(col,col),col);
        );

    );
    return matrix A
)

diagonalizeGW = method()

diagonalizeGW (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    A := beta.matrix;
    D := diagonalize(A);
    return gwClass(D)
    )
