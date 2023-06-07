--This function aims to find the radical of a quadratic space.
load "./diagonalize.m2"


truncateRadical=method()
truncateRadical(Matrix):=(GrothendieckWittClass)=> (A) -> (
    truncatedMatrix= mutableMatrix A;
    if not numRows(A)==numColumns(A) then (
        print ("Input is not a square matrix");
    )
    else (
        truncatedMatrix=mutableMatrix diagonalize(A);
        foundRadical=false;
        for i from 0 to (numRows(A)-1) do (
            if truncatedMatrix_(i, i)==0 then (
                foundRadical=true;
                break
            );
            print("The quadratic space does not have a radical!");
        );
            if foundRadical==true then (
                n=numRows(A)-1;
                for i from 0 to n do (
                B=Matrix(truncatedMatrix);
                truncatedMatrix=mutableMatrix submatrix'(B, {i}, {i});
                    if (n>0) then (n=n-1;)
                    else (break);
            );
            return gwClass(Matrix(truncatedMatrix));
        );
    );
);

A=matrix{{0/1, 0/1, 0/1}, {0/1, 1/1, 0/1}, {0/1, 0/1, 1/1}};
print truncateRadical(A);