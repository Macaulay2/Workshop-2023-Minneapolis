load "./_init.m2"

------------
-- Types
------------

-- We define GrothendieckWittClass to be a new type,
--    meant to represent the isomorphism class of a symmetric bilinear form
--    over a base field.
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

-- This allows us to extract the matrix from a class
matrix GrothendieckWittClass := Matrix => beta -> beta.matrix

-- Check if a constructed class is well-defined
isWellDefined GrothendieckWittClass := Boolean => beta -> (
    
    -- Returns false if a matrix isn't symmetric
    --	  Note that this will also return false if a matrix isn't square,
    --	      so we don't need another check for that.
    if not (transpose(beta.matrix) == beta.matrix) then(
	<< "-- Defining matrix is not symmetric" << endl;
	return false
	);
    
    -- Returns false if a matrix is defined over a field of characteristic two
    if char ring beta.matrix == 2 then(
	<< "-- Package does not support base fields of characteristic two" <<endl;
	return false
	);
    
    -- Returns false if the matrix has determinant zero
    if det beta.matrix == 0 then(
	<< "-- Matrix defines a degenerate bilinear form" <<endl;
	return false
        );
    
    -- Returns false if the defining matrix isn't defined over a field
    if not isField ring beta.matrix then(
	<< "-- Matrix is not defined over a field" <<endl;
	return false
	);
    
    true);

-- Method for taking a direct sum of two Grothendieck-Witt classes
gwAdd = method()

gwAdd(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (
    
    -- Returns an error if the underlying fields of the two inputted symmetric bilinear forms are different
    if not ring beta.matrix === ring gamma.matrix then error "Error: these classes have different underlying fields";

    b := beta.matrix;
    g := gamma.matrix;
    
    return gwClass(b++g)
    )

-- Method for taking a tensor product of two Grothendieck-Witt classes
gwMultiply = method()

gwMultiply(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (

    -- Returns an error if the underlying fields of the two inputted symmetric bilinear forms are different
    if not ring beta.matrix === ring gamma.matrix then error "Error: these classes have different underlying fields";

    b := beta.matrix;
    g := gamma.matrix;
    
    return gwClass(b**g)
    )



-- Method for diagonalizing a square matrix
diagonalize = method()

diagonalize (MutableMatrix) := (MutableMatrix) => (A) -> (
    
    -- TODO return error if not square
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

