

doc ///
    Key
        schubertCodim
        (schubertCodim, Matrix)
        (schubertCodim, List)
    Headline
        Computes the codimension of a schubert determinantal ideal
    Usage
        schubertCodim(M)
        schubertCodim(w)
    Inputs
        w:List
	        or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, outputs the codimension of the corresponding schubert determinantal ideal.
        Example
            schubertCodim(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
            schubertCodim({1,3,2})
///


doc ///
    Key
        matrixSchubertReg
	(matrixSchubertReg, Matrix)
        (matrixSchubertReg, List)
    Headline
        Computes the Castelnuovo-Mumford regularity of a Schubert determinantal ideal
    Usage
    	matrixSchubertReg(M)
        matrixSchubertReg(w)
    Inputs
        w:List
	M:Matrix
    Description
        Text
            Given a partial alternating sign matrix computes the Castelnuovo-Mumford regularity of the corresponding schubert Schubert determinantal ideal by
	    computing the Castelnuovo-Mumford of the antidiagonal initial ideal. Given a permutation in 1-line notation, computes the Castelnuovo-Mumford regularity
	     of the corresponding Schubert determinantal ideal by implementing Theorem 1.1 of CASTELNUOVO–MUMFORD REGULARITY OF MATRIX SCHUBERT VARIETIES by OLIVER PECHENIK, 
	     DAVID E SPEYER, AND ANNA WEIGANDT found at ...
        Example
	    w = {2,3,5,1,4}
	    matrixSchubertReg(w)
	
	    A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}};
            matrixSchubertReg(A)
///



