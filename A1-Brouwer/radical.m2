--This function aims to find the radical of a quadratic space.
GrothendieckWittClass = new Type of HashTable
GrothendieckWittClass.synonym = "Grothendieck Witt Class"

gwClass = method()

-- A class in GW can be constructed from a representing matrix
gwClass (Matrix) := GrothendieckWittClass => M -> (
  new GrothendieckWittClass from {
      symbol matrix => M,
      symbol cache => new CacheTable
      }
)

diagonalize = method()

diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    A := mutableMatrix AnonMut;
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
                rowAdd(A,row,-temp/A_(col,col),col); --more row reduction make every entry below A_(col,col) is zero
                columnAdd(A,row,-temp/A_(col,col),col); --column reduction to keep reduced matrix congruent
            );
        );

    );
    return matrix A 
)



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
            if foundRadical===true then (
                n=numRows(A)-1;
                for i from 0 to n do (
                truncatedMatrix=mutableMatrix submatrix'(matrix truncatedMatrix, {i}, {i});
                    if (n>0) then (n=n-1;)
                    else (break);
            );
            B=matrix truncatedMatrix;
            return gwClass(B);
        );
    );
);
