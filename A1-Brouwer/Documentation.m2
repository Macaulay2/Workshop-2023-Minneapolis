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
            -Returns the Witt Index of q and the matrix representing q_a on V_a as \pm I, after changing basis.
*}





            ---------------------------- documentation starts here-------------

beginDocumentation()


undocumented {
-- put undocumented commands here.
    }

document {
	Key => GrothendieckWitt,
	Headline => "Package for calculating Grothendieck-Witt groups and calculating local and global A1-degrees", 
	Acknowledgement => {}
}

document {
	Key => {(diagonalize, Matrix), diagonalize},
	Headline => "diagonalizing a symmetric matrix via congruence",
	Usage => "diagonalize(M)",
	Inputs => {
		Matrix => "M" => {"a symmetric matrix over any field"}
		},
	Outputs => { Matrix => { "a diagonal matrix congruent to", TT "M" }},
	PARA {"Given a symmetric matrix ", TT "M", " over any field, this command gives a diagonal matrix congruent to", TT "M",".",
	"Note that the order in which the diagonal terms appear is not specified."},
	EXAMPLE lines ///
		 M=matrix(ZZ/17, {{7, 9}, {9, 6}});
		 result=diagonalize(M);
		 assert(result===matrix(ZZ/17, {{7, 0}, {0, -8}}));
	 	 ///,
     	}

document {
	Key => {(diagonalForm, GrothendieckWittClass), GrothendieckWittClass},
	Headline => "nice diagonal representative of a Grothendieck-Witt Class",
	Usage => "diagonalize(G)",
	Inputs => {
		GrothendieckWittClass => "G" => {"a Grothendieck-Witt Class over", TEX ///$\mathbb{R}$///, "or",
		TEX///$\mathbb{C}$///, "."},
	},
	Outputs => {GrothendieckWittClass => {"a Grothendieck-Witt Class isomorphic to", TT "G",
	"whose matrix representative is diagonal."}},
	PARA {"Given a Grothendieck-Witt Class," TT "G", "over", TEX///$\mathbb{R}$///, "or", TEX///$\mathbb{C}$///,
	"returns a Grothendieck-Witt Class whose representing matrix has the following form: the first several diagonal blocks",
	"are of the form" TEX///$\begin{bmatrix}1 & 0 \\ 0 & 1\end{bmatrix}$///, ", and the number of such blocks equals the number of",
	"hyperbolic planes contained in the quadratic space; ",
	"the rest forms an identity matrix whose rank equals the dimension of the anisotropic part of the quadratic space."},
	EXAMPLE lines ///
	M=matrix(RR, {{1, 2, 3, 4}, {2, 4, 5, 16}, {3, 5, 7, 8}, {4, 16, 8, 19}});
	G=gwClass(M);
	M'=diagonalForm(G);
	A=matrix(RR, {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}});
	assert(M6.matrix===A);
	///,
	EXAMPLE lines ///
	M=matrix(CC, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
	G=gwClass(M);
	M'=diagonalForm(G);
	assert(M'.matrix===matrix(CC, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
	///,
}

document {
    Key => {(diagonalizeOverInt, Matrix), diagonalizeOverInt},
	Headline => "diagonalizing a symmetric matrix with integer entries via congruence",
	Usage => "diagonalizeOverInt(M)",
	Inputs => {
		Matrix => "M" => {"a symmetric matrix with integer entries"}
		},
	Outputs => { Matrix => { "a diagonal matrix congruent to ", TT "M" }},
	PARA {"Given a symmetric matrix ", TT "M", " over the integers, this command gives a diagonal matrix congruent to", TT "M","."},
	EXAMPLE lines ///
		M=matrix(ZZ, {{6, 5},{5, 9}});
		A=matrix(ZZ, {{6, 0}, {0, 174}});
		assert(diagonalizeOverInt(testMatrix6)===A);
	 	 ///,
        }

document {
    Key => {(diagonalizeCB, Matrix), diagonalizeCB},
	Headline => "Row/column reducing a square matrix to a diagonal form within an n x k matrix",
	Usage => "diagonalizeCB(M)",
	Inputs => {
		Matrix => "M" => {"an n x k matrix"}
		},
	Outputs => { Matrix => { "the n x k matrix whose first n columns have been row/column reduced to a diagonal matrix" }},
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
	Headline => "Witt decomposition of an invertible matrix not over the reals or complex numbers",
	Usage => "wittDecomp(M)",
	Inputs => {
		Matrix => "M" => {"an invertible matrix over any field that isn't the real or complex numbers"}
		},
	Outputs => { (ZZ,Matrix) => { "the number of hyperbolic forms that make up the bilinear form represented by ", TT "M", ", and the anisotropic part of the form represented as a matrix" }},
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
	Headline => "Grothendieck Witt class of a symmetric matrix",
	Usage => "gwClass(M)",
	Inputs => {
		Matrix => "M" => {"a symmetric matrix"}
		},
	Outputs => { GrothendieckWittClass => { "the isomorphism class of a symmetric bilinear form represented by ", TT "M" }},
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
	Headline => "base field of a Grothendieck Witt class",
	Usage => "baseField(beta)",
	Inputs => {
		GrothendieckWittClass => "beta" => {"the isomorphism class of a symmetric bilinear form"}
		},
	Outputs => { Ring => { "the base field of the Grothendieck-Witt class ", TT "beta" }},
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
    Key => {(squarefreePart, QQ), (squarefreePart, ZZ), squarefreePart},
	Headline => "smallest magnitude representative of a square class over the rationals or integers"
	Usage => "squarefreePart(q)
				squarefreePart(n)",
	Inputs => {
		QQ => "q" => {a rational number},
		ZZ => "n" => {an integer}
		},
	Outputs => { ZZ => { "the smallest magnitude integer in the square class of ", TT "n"}},
	PARA {"Given a rational number (or integer), ", TT "q", ", this command outputs the smallest magnitude integer, ",
                TEX///$m$///, ", such that ", TEX///$q=lm$///, " for some rational number (or integer) ",
				TEX///$l$///, "."},
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
    Key => {(localAlgebraBasis, List, Ideal), localAlgebraBasis},
	Headline => "basis for local k-alegbra ", TEX///$Q_p(f)$///,
	Usage => "localAlgebraBasis(L,p)",
	Inputs => {
		List => "L" => {"list of polynomials ", TEX///$f=(f_1, \dots ,f_n)$///, " over the same ring"},
		Ideal => "p" => {"prime ideal of an isolated zero"}
	},
	Outputs => { List => {"a list of basis elements of the local k-algebra", TEX///$Q_p(f)$/// }},
	PARA {"Given an endomorphism of affine space, ", TEX///$f=(f_1,\dots ,f_n)$///,
			", given as a list of polynomials called ", TT "L", " and the prime ideal of an isolated zero, this command returns a list of basis elements of the local k-algebra", TEX///$Q_p(f)$///, "."}
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
	Inputs => {"An endomorphism", TEX///$f=(f_{1},...,f_{n}):\mathbb{A}^{n}_{k}\to\mathbb{A}^{n}_{k}$///, "as a list", TT "f={f_1, ..., f_n}"},
	Outputs => { gwClass => { "the Grothendieck-Witt class defining the", TEX///$\mathbb{A}^{1}$///, "-degree of the endomorphism" TT "f", "." }},
	PARA {"Given an endomorphism of affine space, ", TEX///$f=(f_{1},...,f_{n}):\mathbb{A}^{n}_{k}\to\mathbb{A}^{n}_{k}$///, "as a list", TT "f={f_1, ..., f_n}", ", this method computes the", TEX///$\mathbb{A}^{1}$///, "-degree via Bezoutians."},
	EXAMPLE lines ///
		 S = QQ[z]
		 f = {z^4 + z^3 - z^2 - z}
		 globalA1Degree(f)
	 	 ///,
        }
    
document {
    Key => {(locallA1Degree, List, Ideal), localA1Degree},
    Headline => "Bezoutian for local endomorphism f",     
	Usage => "localA1Degree(L, p)",
	Inputs => {
	    List => "L" => {"List of polynomials ", TEX///$f=(f_{1},...,f_{n})$///, " over ring k defining an endomorphism ", 
		TEX///$f:\mathbb{A}^{n}_{k}\to\mathbb{A}^{n}_{k}$///}
	    Ideal => "p" => {"prime ideal of an isolated zero for " TT "f"}
	    }, 
	Outputs => { 
	    Matrix => { "the Bezoutian defining the local ", TEX///$\mathbb{A}^{1}$///, "-degree of the endomorphism" TT "f", 
		" for local ring R_(p), where R is the ring over which", TT "f", " is defined." }},
	PARA {"For a local ring " TT k$ ", given an endomorphism of affine space, ", TEX///$f=(f_{1},...,f_{n}):\mathbb{A}^{n}_{k}\to\mathbb{A}^{n}_{k}$///, "as a list", TT "f={f_1, ..., f_n}", ", this method computes the local ", TEX///$\mathbb{A}^{1}$///, "-degree via Bezoutians."},
	EXAMPLE lines ///
		 S = QQ[z]
		 p = ideal (z) 
		 L = {z^4 + z^3 - z^2 - z}
		 localA1Degree(L, p)
	 	 ///,
		 lines ///
		 S=QQ[x,y]
		 L = {x^2+y^2-1, x^2-y-1}
		 p = ideal (x, y)
		 localA1Degree(L, p)
	 	 ///,
        }
