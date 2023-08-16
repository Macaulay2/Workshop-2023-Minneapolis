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
    k := ring A;
    if not (k === RR or instance(kk,RealField) or k === QQ) then(
        error "Only implemented over QQ and RR";
        );
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
    k := ring A;
    if not (k === RR or instance(kk,RealField) or k === QQ) then(
        error "Only implemented over QQ and RR";
        );
    negDiagEntries := 0;
    for i from 0 to (numRows(A)-1) do (
        if A_(i,i) < 0 then (
            negDiagEntries = negDiagEntries + 1;
            );
        );
    return(negDiagEntries);
    )
