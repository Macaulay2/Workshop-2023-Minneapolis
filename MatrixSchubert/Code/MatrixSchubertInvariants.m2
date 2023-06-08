------------------------------------------------------------
--code for exploring invariants of MatrixSchubert varieties
------------------------------------------------------------


------------------------------------------
--INPUT: matrixSchubertReg, takes a permutation or an ASM 
--       and an optional strategy computation argument. The
--       strategy options are ADI (computes the antiDiagonalInitial
--       ideal and then computes the Castelnuovo-Mumford regularity)
--       and PSW (computes the Castelnuovo-Mumford regularity of the
--       matrix Schubert variety using Theorem 1.1 of Pechenki-Speyer-Weigandt)
--OUTPUT: returns the Castelnuovo-Mumford reguarity of the matrix 
--        Schubert variety
------------------------------------------

matrixSchubertRegADI = method()
matrixSchubertRegADI List := ZZ => (w) -> (
    if not (isPerm w) then error ("Expecting a permutation.");
    
    I := antiDiagInit w;
    if I == 0 then return 0;
    return regularity(I) -1;
    
);

matrixSchubertRegWPS = method()
matrixSchubertRegWPS List := ZZ => (w) -> (

     if not (isPerm w) then error ("Expecting a permutation.");
     
     return rajIndex(w) - permLength(w);
         
);

-*
matrixSchubertReg() = method(
    Options => {
	Strategy => ADI
	}    
)

matrixSchubertReg (List) := opts -> 
*-

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation or a partial ASM
--OUTPUT: the Castlenuovo-Mumford regularity of I_A or I_w
----------------------------------------
schubertReg= method()
schubertReg Matrix := ZZ => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    regularity(antiDiagInit A)
    );
schubertReg List := ZZ => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
     regularity(antiDiagInit w)
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
