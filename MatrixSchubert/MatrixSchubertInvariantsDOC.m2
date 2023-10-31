
doc ///
    Key
        schubReg
        (schubReg, List)
	(schubReg, Matrix)
    Headline
        compute the Castelnuovo-Mumford regularity of the quotient by a Schubert determinantal ideal or ASM ideal
    Usage
        schubReg(w)
        schubReg(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, computes the Castelnuovo-Mumford regularity of the quotient by the corresponding alternating sign matrix ideal or Schubert determinantal ideal.  In the case of a permutation in 1-line notation, computes the Castelnuovo-Mumford regularity of the corresponding Schubert determinantal ideal by implementing Theorem 1.1 of 

            @UL {
                {"Oliver Pechenik, David Speyer, and Anna Weigandt, ", EM "Castelnuovo-Mumford regularity of matrix Schubert varieties, ", arXiv  "2111.10681"}
            }@
	
	    In the case of a partial permutation, computes the regularity using the antidiagonal initial ideal, a valid strategy in light of 
	    
	     @UL {
                {"Aldo Conca and Matteo Varbaro, ", EM "Square-free Gröbner degenerations, ", arXiv  "1805.11923", ", ", "Invent. Math. 221 (2020), no. 3, 713–730."}
            }@
	    
        Example
            w = {2,3,5,1,4}
            schubReg(w)
            A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}};
            schubReg(A)
///


doc ///
    Key
        schubCodim
        (schubCodim, List)
        (schubCodim, Matrix)
    Headline
        compute the codimension (i.e., height) of a Schubert determinantal ideal or ASM ideal
    Usage
        schubCodim(w)
        schubCodim(M)
    Inputs
        w:List
	        or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, outputs the codimension of the corresponding alternating sign matrix ideal or Schubert determinantal ideal.
        Example
            schubCodim(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
            schubCodim({1,3,2})
///


doc ///
    Key
    	isSchubCM
    	(isSchubCM, Matrix)
    	(isSchubCM, List)
    Headline
    	whether an ASM variety is Cohen-Macaulay
    Usage
    	isSchubCM(A)
    	isSchubCM(w)
    Inputs
	A:Matrix
		or {\tt w} is a @TO List@
    Outputs
	:Boolean
    Description
    	Text
	    Given an alternating sign matrix $A$ (resp. permutation $w$),
	    checks whether $R/I_A$ (resp. $R/I_w$) is Cohen-Macaulay.
	    If the input is a permutation $w$, the output is always true
	    since $I_w$ is a Schubert determinantal ideal, and a theorem
	    of Fulton says $R/I_w$ is always Cohen-Macaulay.
        Example
            A = matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}}
	    isSchubCM(A)
	    w = {1,3,2}
            isSchubCM(w)
///

doc ///
    Key 
        KPolynomialASM
        (KPolynomialASM, Matrix)
    Headline
        compute the K Polynomial of an ASM variety
    Usage
        KPolynomialASM(A)
    Inputs 
        A:Matrix
    -- Outputs
        -- :PolynomialRing --Not sure what to put here, PolynomialRing makes the compiler angry
    Description
        Text
            Given a partial ASM $A$, compute the K Polynomial of its corresponding Ideal, defined as the numerator of its Hilbert series. The multi-degree variables are indexed along rows.
        Example
            A = matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}};
            KPolynomialASM(A)
///


-- doc ///
--     Key
--         (matrixSchubertRegADI, List)
--         matrixSchubertRegADI
--     Headline
--         compute the Castelnuovo-Mumford regularity of a Schubert determinantal ideal using antidiagonal initial ideal
--     Usage
--         matrixSchubertRegADI(w)
--     Inputs
--         w: List 
--     Description
--         Text
--             Given a permutation in 1-line notation, computes the regularity of the corresponding Schubert determinantal ideal by computing it for the antidiagonal initial ideal.
--         Example 
--             w = {2,3,5,1,4}
--             matrixSchubertReg(w)
-- ///