------------------------------------------------------------
--code for exploring invariants of MatrixSchubert varieties
------------------------------------------------------------

------------------------------------------
--INPUT: matrixSchubertRegADI, takes a permutation in 1-line notation
--OUTPUT: returns the Castelnuovo-Mumford reguarity of the matrix 
--        Schubert variety by computing the regularity of the antidiagonal initial ideal
------------------------------------------
-- matrixSchubertRegADI = method()
-- matrixSchubertRegADI List := ZZ => (w) -> (
--     if not (isPerm w) then error ("Expecting a permutation.");
    
--     I := antiDiagInit w;
--     if I == 0 then return 0;
--     return regularity(I) -1;   
-- )

matrixSchubertReg = method()
matrixSchubertReg List := ZZ => (w) -> (
     if not (isPerm w) then error ("Expecting a permutation.");
     return rajIndex(w) - permLength(w);
)
matrixSchubertReg Matrix := ZZ => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    -- TODO: Check if Matrix is permutation matrix
        -- if it is, use rajindex formula
    I := antiDiagInit A;
    if I == 0 then return 0;
    return regularity(I) -1;
);

schubertCodim = method() 
schubertCodim Matrix := ZZ => A -> (
    if not (isPartialASM A) then error("The input must be a partial alternating sign matrix");
    codim antiDiagInit A
)
schubertCodim List := ZZ => (w) -> (
    if not (isPerm w) then error("The input must be a permutation in one line notation");
    permLength w
)