-Type: GrothendieckWittClass

    (1). Matrix -> GrothendieckWittClass
                gwClass (Matrix)
         Matrix <- GrothendieckWittClass 
                (GrothendieckWittClass).matrix
    (2). Methods:
        (a). baseField (GrothendieckWittClass) -> (Ring)
                Returns base field of a GrothendieckWittClass
        (b). diagonalize (Matrix) -> (Matrix)
                Input a matrix, returns the diagonalized form of the matrix.
        (c). gwAdd (GrothendieckWittClass, GrothendieckWittClass) -> GrothendieckWittClass
                Returns the direct sum of two GrothendieckWittClasses (isomorphism class of direct sum of two quadratic spaces)
        (d). gwMultiply (GrothendieckWittClass, GrothendieckWittClass) -> GrothendieckWittClass
                Returns the tensor product of two GrothendieckWittClasses (isomorphism class of direct sum of two quadratic spaces)
        (e). isWellDefined (GrothendieckWittClass) -> (Boolean)
                Returns true if matrix representing a GW class is defined over a field k, char(k)\neq 2, symmetric, nonsingular
-Functions:
    (1). Basic Matrix Properties Check
        (a). isSquare (Matrix)-> (Boolean)
        (b). isSquareAndSymmetric (Matrix) -> (Boolean)
        (c). isUpperLeftTriangular (Matrix) -> (Boolean)
        (d). isDiagonal  (Matrix) -> (Boolean)
    (2). easyIsomorphicGW (GrothendieckWittClass, GrothendieckWittClass) -> Boolean
         easyIsomorphicGW (GrothendieckWittClass, GrothendieckWittClass, HeightBound=c) -> Boolean
            Input two matrices (over the same base field), returns true if they are congruent by matrix whose entries have height up to HeightBound.
            *** c takes values in ZZ
            *** Default HeightBound is 3
    (3). easyupperlefttriangular (Matrix) -> (Boolean)
            Returns true if it is upper left triangular after permutation of columns.
    (4). rationalsimplify (Matrix) -> (Matrix)
            Diagonalizes a matrix, making entries square-free, returns matrix after splitting off hyperbolic part
    (5). splitOffObviousHyperbolic (Matrix) -> (ZZ, Matrix)
            Input representing matrix of quadratic space (V, q);
            Returns 1 if quadratic space (V, q) contains hyperbolic plane(s), 0 otherwise;
            Returns complementing quadratic space to the hyperbolic plane(s)
    (6). squarefreepart (QQ) -> (ZZ)
         squarefreepart (ZZ) -> (ZZ)
            Input rational number, returns maximal positive square-free factor of product of denominator and numerator.
            Input integer, returns maximal positive square-free factor.
    (7). wittDecomp (Matrix, Ring (field) ) -> (ZZ, Matrix) 
            -(V, q) quadratic space over a field F, quadratic form represented by matrix M.
            -Witt decomposition: (V_t, q_t)\perp (V_h, q_h)\perp (V_a, q_a)
                (respectively: totally isotropic, hyperbolic, anisotropic)
            -Returns the Witt Index of q and the matrix representing q_a on V_a.
    (8). wittDecompInexact (Matrix, InexactFieldFamily) -> (ZZ, Matrix)
            -InexactFieldFamily must be RR or CC
            -Returns the Witt Index of q and the matrix representing q_a on V_a as \pm I, after changing basis.





            ---------------------------- documentation starts here-------------

beginDocumentation()

document {
	Key => {(diagonalize, Matrix), diagonalize},
	Usage => "diagonalize(M)",
	Inputs => {"M"},
	Outputs => { Matrix => { "a diagonal matrix congruent to", TT "M" }},
	PARA {"Given a symmetric matrix ", TT "M", " over any field, this command gives a diagonal matrix congruent to", TT "M","."},
	EXAMPLE lines ///
		 R = QQ[x,y]
		 F = {y^2-x^2-1,x-y^2+4*y-2}
		 I = ideal F
		 regularRep(y,I)
		 S = R/I
		 regularRep(y)
	 	 ///,
     	}

document {
        Key => {(diagonalizeOverInt, Matrix), diagonalizeOverInt},
	Usage => "diagonalizeOverInt(M)",
	Inputs => {"M"},
	Outputs => { Matrix => { "a diagonal matrix congruent to ", TT "M", "." }},
	PARA {"Given a symmetric matrix ", TT "M", " over the integers, this command gives a diagonal matrix congruent to", TT "M","."},
	EXAMPLE lines ///
		 R = QQ[x,y]
		 F = {y^2-x^2-1,x-y^2+4*y-2}
		 I = ideal F
		 regularRep(y,I)
		 S = R/I
		 regularRep(y)
	 	 ///,
        }

document {
        Key => {(wittDecomp, Matrix), wittDecomp},
	Usage => "wittDecomp(M)",
	Inputs => {"M"},
	Outputs => { (ZZ,Matrix) => { "the number of hyperbolic forms that make up the bilinear form represented by ", TT "M", ", and the anisotropic part of the form represented as a matrix." }},
	PARA {"Given an invertible matrix ", TT "M", ", which represents a non-degenerate bilinear form, ", TT "b"," , this command gives us ", TT "(n,A)", " where ", TT "b="},
	EXAMPLE lines ///
		 R = QQ[x,y]
		 F = {y^2-x^2-1,x-y^2+4*y-2}
		 I = ideal F
		 regularRep(y,I)
		 S = R/I
		 regularRep(y)
	 	 ///,
        }