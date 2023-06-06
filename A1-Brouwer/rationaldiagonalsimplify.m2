--takes in diagonal matrix over QQ, removes squares from entries and splits off obvious hyperbolics

load "matrixBooleans.m2"
load "squarefreepart.m2"
load "splitoffobvioushyperbolics.m2"

rationalDiagonalSimplify = method()

rationalDiagonalSimplify (Matrix) := (ZZ,Matrix) => (A) -> (
    B = mutableMatrix A;
    for i from 0 to (numRows(B)-1) do (
        B_(i,i) = squarefreePart(B_(i,i));
        );
    C = matrix B;
    print(B);
    return(splitOffObviousHyperbolics(C))
    )