

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
