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

schubReg = method()
schubReg List := ZZ => w -> (
     if not (isPerm w) then error ("The input must be a partial alternating sign matrix or a permutation.");
     return rajIndex(w) - permLength(w);
)
schubReg Matrix := ZZ => A -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    -- TODO: Check if Matrix is permutation matrix
        -- if it is, use rajindex formula
    I := antiDiagInit A;
    if I == 0 then return 0;
    return regularity(I) -1;
);

schubCodim = method() 
schubCodim Matrix := ZZ => A -> (
    if not (isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    codim antiDiagInit A
)
schubCodim List := ZZ => w -> (
    if not (isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    permLength w
)


--------------------------------------------------
--**KPolynomialASM**-
--input: a partial ASM A
--output: the Kpolynomial of ASM variety for A
     --using multidegree where variables are indexed along rows
--TODO: add option for doubly graded
--------------------------------------------------
KPolynomialASM = method()
KPolynomialASM Matrix := ZZ => A -> (
    I := schubertDetIdeal(A,CoefficientRing=>ZZ/2);
    R := ring I;
    kk := coefficientRing R;
    possibleDegs := apply(numrows A, i-> toList insert(i,1,(numrows(A)-1):0));
    degs := splice apply(possibleDegs, i->(numrows(A):i));
    Q := kk[R_*, Degrees => degs];
    numerator hilbertSeries sub(I,Q)
);
