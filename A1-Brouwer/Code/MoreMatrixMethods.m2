numNonzeroDiagEntries = method()
numNonzeroDiagEntries (Matrix) := (Matrix) => (A) -> (
    nonzeroDiagEntries := 0;
    for i from 0 to (numRows(A)-1) do (
        if A_(i,i) != 0 then (
            nonzeroDiagEntries = nonzeroDiagEntries + 1;
            );
        );
    return(nonzeroDiagEntries);
    )

numPosDiagEntries = method()
numPosDiagEntries (Matrix) := (Matrix) => (A) -> (
    posDiagEntries := 0;
    for i from 0 to (numRows(A)-1) do (
        if A_(i,i) > 0 then (
            posDiagEntries = posDiagEntries + 1;
            );
        );
    return(posDiagEntries);
    )

numNegDiagEntries = method()
numNegDiagEntries (Matrix) := (Matrix) => (A) -> (
    negDiagEntries := 0;
    for i from 0 to (numRows(A)-1) do (
        if A_(i,i) < 0 then (
            negDiagEntries = negDiagEntries + 1;
            );
        );
    return(negDiagEntries);
    )
