--This function aims to find the radical of a quadratic space.
load "diagonalize.m2"
load "GW-type.m2"




-- This is reliant on the diagonalize() method being applicable for singular matrices

truncateRadical=method()
truncateRadical(Matrix):=(Matrix)=> (A) -> (
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
            if foundRadical===true then (
                n=numRows(A)-1;
                for i from 0 to n do (
                truncatedMatrix=mutableMatrix submatrix'(matrix truncatedMatrix, {i}, {i});
                    if (n>0) then (n=n-1;)
                    else (break);
            );
            B=matrix truncatedMatrix;
            return B;
        );
    );
);
