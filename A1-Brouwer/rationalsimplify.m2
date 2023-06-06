--takes in symmetric matrix over QQ and diagonalizes, removes squares from entries, and splits off hyperbolic forms that immediately appear as <a> + <-a>

load "diagonalize.m2"
load "matrixBooleans.m2"
load "squarefreepart.m2"
load "splitoffobvioushyperbolics.m2"

rationalSimplify = method()

rationalSimplify (Matrix) := (ZZ,Matrix) => (A) -> (
    --matrix must be square
    if numRows(A) != numColumns(A) then (
	error "Matrix is not square.";
	);
    --matrix must be nonsingular
    if det(A) == 0 then (
	error "Matrix is singular.";
	);
    --matrix must be symmetric
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    --diagonalize the matrix
    B = mutableMatrix(diagonalize(A));
    --replace entry with smallest magnitude integer in square class
    for i from 0 to (numRows(B)-1) do (
        B_(i,i) = squarefreePart(B_(i,i));
        );
    C = matrix B;
    --split off hyperbolic forms <a> + <-a>
    return(splitOffObviousHyperbolics(C))
    )