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