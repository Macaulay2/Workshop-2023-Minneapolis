--Tom Hagedorn, Joel Louwsma

--given diagonal matrix, split off any <a>+<-a> and return number of times we can do this as well as smaller matrix with none of these

splitOffObviousHyperbolic = method()

splitOffObviousHyperbolic (Matrix) := (ZZ,Matrix) => (A) -> (
    --matrix must be symmetric
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    foundHyperbolic = 0;
    remainingMatrix = A;
    for i from 0 to (numRows(A)-1) do (
        for j from (i+1) to (numRows(A)-1) do (
            if (A_(i,i) == -A_(j,j) and A_(i,i) != 0) then (
                 foundHyperbolic = 1;
                 remainingMatrix = submatrix'(A,{i,j},{i,j});
                 return(foundHyperbolic,remainingMatrix)
                 );
            );
        );
    return(foundHyperbolic,remainingMatrix);
    )


splitOffObviousHyperbolics = method()

splitOffObviousHyperbolics (Matrix) := (ZZ,Matrix) => (A) -> (
    numberHyperbolics=0;
    notFinished=1;
    while (notFinished == 1) do (
        currentState = splitOffObviousHyperbolic(A);
        notFinished = currentState_0;
        numberHyperbolics = numberHyperbolics + notFinished;
        A = currentState_1;
        );
    return (numberHyperbolics,A);
    )
