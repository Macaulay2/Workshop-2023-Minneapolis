
doc ///
    Key
        (matrixSchubertReg, List)
	    (matrixSchubertReg, Matrix)
        matrixSchubertReg
    Headline
        compute the Castelnuovo-Mumford regularity of a Schubert determinantal ideal
    Usage
        matrixSchubertReg(w)
        matrixSchubertReg(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a partial alternating sign matrix computes the Castelnuovo-Mumford regularity of the corresponding schubert Schubert determinantal ideal by computing the Castelnuovo-Mumford of the antidiagonal initial ideal. Given a permutation in 1-line notation, computes the Castelnuovo-Mumford regularity of the corresponding Schubert determinantal ideal by implementing Theorem 1.1 of 

            @UL {
                {"Oliver Pechenik, David Speyer, and Anna Weigandt, ", EM "Castelnuovo-Mumford regularity of matrix Schubert varieties, ", arXiv  "2111.10681"}
            }@
        Example
            w = {2,3,5,1,4}
            matrixSchubertReg(w)
            A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}};
            matrixSchubertReg(A)
///


doc ///
    Key
        (matrixSchubertCodim, List)
        (matrixSchubertCodim, Matrix)
        matrixSchubertCodim
    Headline
        compute the codimension of a schubert determinantal ideal
    Usage
        matrixSchubertCodim(w)
        matrixSchubertCodim(M)
    Inputs
        w:List
	        or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, outputs the codimension of the corresponding schubert determinantal ideal.
        Example
            matrixSchubertCodim(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
            matrixSchubertCodim({1,3,2})
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
