diagonalize = method()

diagonalize (MutableMatrix) := (MutableMatrix) => (A) -> (
    n=numRows(A);
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
            if A_(col,col)==0 then (print "Error: Matrix A was singular"; return A;);
        );
        --entry in A_(col,col) is non-zero at this point
         for row from (col+1) to (n-1) do (
            temp=A_(row,col);
            rowAdd(A,row,-temp/A_(col,col),col);
            columnAdd(A,row,-temp/A_(col,col),col);
        );

    );
    return A
)

A=matrix{{1/1, 0},{0 ,1}};
Amut= mutableMatrix A;
B=matrix{{0,1/1},{1 ,0 }};
Bmut= mutableMatrix B;
C=matrix{{1,1/1},{1 ,1 }};
Cmut= mutableMatrix C;
D=matrix{{-1/1,-1,1,1},{-1,1,1,0 },{1,1,0,0},{1,0,0,0}};
Dmut= mutableMatrix D;

print diagonalize(Amut)
print diagonalize(Bmut)
print diagonalize(Cmut)
print diagonalize(Dmut)


