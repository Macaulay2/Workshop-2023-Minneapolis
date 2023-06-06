needsPackage "RationalPoints2"
--diagonalize method
--given a symmetric invertible matrix, this function outputs a diagonal matrix congruent to og matrix
diagonalize = method()

diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    A := mutableMatrix AnonMut;
    n:=numRows(A);
    for col from 0 to (n-1) do (
        if A_(col,col) == 0 then ( --if diagonal entry in column "col" is zero
            for row from col+1 to n-1 do ( 
                if A_(row,col) != 0 then ( --scan for nonzero entries below the diagonal entry
                    rowAdd(A,col,1,row); --
                    columnAdd(A,col,1,row);
                    break;
                );
            );
            if A_(col,col)==0 then (error "Error: Matrix A was singular"; return A;);
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

A=matrix{{1/1, 0},{0 ,1}};
B=matrix{{0,1/1},{1 ,0 }};
C=matrix{{1,1/1},{1 ,1 }};
D=matrix{{-1/1,-1,1,1},{-1,1,1,0},{1,1,0,0},{1,0,0,0}};

print diagonalize(A)
print diagonalize(B)
print diagonalize(C)
print diagonalize(D)


