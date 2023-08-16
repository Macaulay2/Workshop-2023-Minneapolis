nondegeneratePartDiagonal = method()
nondegeneratePartDiagonal (Matrix) := (Matrix) => (A) -> (
    diagA := congruenceDiagonalize(A);
    i := 0;
    while (i < numRows(diagA)) do (
        if (diagA_(i,i) == 0) then (
            diagA = submatrix'(diagA,{i},{i});
            )
        else (
            i = i+1;
            );
        );
    return (diagA);
    )

-- Boolean returning whether a symmetric bilinear form or Grothendieck-Witt class is nondegenerate	
isNondegenerate = method()
isNondegenerate (Matrix) := (Boolean) => (A) -> (
    return (numRows(A) == numRows(nondegeneratePartDiagonal(A)));
    )

isNondegenerate (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    return (isNondegenerate(alpha.matrix));
    )

-- Boolean returning whether a symmetric bilinear form or Grothendieck-Witt class is degenerate	
isDegenerate = method()
isDegenerate (Matrix) := (Boolean) => (A) -> (
    return (not isNondegenerate(A));
    )

isDegenerate (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    return (not isNondegenerate(alpha));
    )
