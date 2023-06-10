{*-Type: GrothendieckWittClass

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
            -Returns the Witt Index of q and the matrix representing q_a on V_a as \pm I, after changing basis.*}





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
        Key => {(diagonalizeCB, Matrix), diagonalizeCB},
	Usage => "diagonalizeCB(M)",
	Inputs => {"M"},
	Outputs => { Matrix => { "the n x k matrix whose first n columns have been row/column reduced to a diagonal matrix." }},
	PARA {"Given an n x k matrix, ", TT "M", ", this command outputs the n x k matrix whose first n columns have been row/column reduced to a diagonal matrix."},
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
	PARA {"Given an invertible matrix, ", TT "M", ", over any field that isn't the real or complex numbers which represents a non-degenerate bilinear form, ",TEX///$\beta$///,
                ", this commmand outputs", TEX///$n$ ///, " and ", TEX///$A$ ///, " where", 
                TEX/// $\beta = n\mathbb{H} + \beta_a///, ", and ", TEX///$A$///, " represents ", TEX///$\beta_a///, "."},
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
        Key => {(gwClass, Matrix), gwClass},
	Usage => "gwClass(M)",
	Inputs => {"M"},
	Outputs => { GrothendieckWittClass => { "the isomorphism class of a symmetric bilinear form represented by ", TT "M", "." }},
	PARA {"Given a symmetric matrix, ", TT "M", ", this command outputs an object of type ", TT "GrothendieckWittClass", ". ",
                "This output has the representing matrix, ", TT "M", ", and the base field of the matrix stored in its CacheTable."},
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
        Key => {(baseField, GrothendieckWittClass), baseField},
	Usage => "baseField(beta)",
	Inputs => {"beta"},
	Outputs => { Ring => { "the base field of the Grothendieck-Witt class ", TT "beta", "." }},
	PARA {"Given the isomorphism class of a symmetric bilinear form, ", TT "beta", 
                ", this command outputs the base field of the form."},
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
	PARA {"Given an invertible matrix, ", TT "M", ", over any field that isn't the real or complex numbers which represents a non-degenerate bilinear form, ",TEX///$\beta$///,
                ", this commmand outputs", TEX///$n$ ///, " and ", TEX///$A$ ///, " where", 
                TEX/// $\beta = n\mathbb{H} + \beta_a///, ", and ", TEX///$A$///, " represents ", TEX///$\beta_a///, "."},
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
        Key => {(squarefreePart, QQ or ZZ), squarefreePart},
	Usage => "squarefreePart(n)",
	Inputs => {"n"},
	Outputs => { ZZ => { "the smallest magnitude integer in the square class of ", TT "n", "." }},
	PARA {"Given a rational number, ", TT "n", ", this command outputs the smallest magnitude integer, ",
                TEX///$m$///, ", such that "},
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
        Key => {(globalA1Degree, List), globalA1Degree},
	Usage => "globalA1Degree(f)",
	Inputs => {"An endomorphism" TEX///$f=(f_{1},...,f_{n}):\mathbb{A}^{n}_{k}\to\mathbb{A}^{n}_{k}$///, "as a list" TT"f={f_1, ..., f_n}"},
	Outputs => { gwClass => { "the Grothendieck-Witt class defining the" TEX///$\mathbb{A}^{1}$/// "degree of the endomorphism" TT "f", "." }},
	PARA {"Given an endomorphism of affine space, ", TT "n", ", this command outputs the smallest magnitude integer, ",
                TEX///$m$///, ", such that "},
	EXAMPLE lines ///
		 R = QQ[x,y]
		 F = {y^2-x^2-1,x-y^2+4*y-2}
		 I = ideal F
		 regularRep(y,I)
		 S = R/I
		 regularRep(y)
	 	 ///,
        }