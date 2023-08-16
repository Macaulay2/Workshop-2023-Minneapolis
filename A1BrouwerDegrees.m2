--A1BrouwerDegrees.m2
newPackage(
    "A1BrouwerDegrees",
    Version=>"1.0",
    Date=>"June 5, 2023",
    Authors=>{
        {Name=>"Nikita Borisov",
	 Email=>"nborisov@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/nikita-borisov"},
        {Name=>"Thomas Brazelton",
	 Email=>"tbraz@math.upenn.edu",
	 HomePage=>"https://www2.math.upenn.edu/~tbraz/"},
        {Name=>"Frenly Espino",
	 Email=>"frenly@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/frenly-espino"},
         {Name=>"Tom Hagedorn",
	 Email=>"hagedorn@tcnj.edu",
	 HomePage=>"https://hagedorn.pages.tcnj.edu/"},
        {Name=>"Zhaobo Han",
	 Email=>"zbtomhan@sas.upenn.edu",
	 HomePage=>"https://www.linkedin.com/in/zhaobo-han-77b1301a2/"},
     	{Name=>"Jordy Lopez Garcia",
	 Email=>"jordy.lopez@tamu.edu",
	 HomePage=>"https://jordylopez27.github.io/"},
        {Name=>"Joel Louwsma",
	 Email=>"jlouwsma@niagara.edu",
	 HomePage=>"https://www.joellouwsma.com/"},
        {Name=>"Andrew Tawfeek",
	 Email=>"atawfeek@uw.edu",
	 HomePage=>"https://www.atawfeek.com/"},
        {Name=>"Wern Juin Gabriel Ong",
	 Email=>"gong@bowdoin.edu",
	 HomePage=>"https://wgabrielong.github.io/"}
	},
    Headline=>"TODO",
    PackageImports=>{
	"Parametrization",
	"RealRoots",
	"RationalPoints2"
	},
    PackageExports=>{},
    DebuggingMode=>true
    )


export{
    -- ArithmeticMethods.m2
    "squarefreePart",
    "legendreBoolean",
    "localAlgebraBasis",
    "primeFactors",
    "PadicValuation",
    "squareSymbol",

    --MatrixMethods.m2
    "congruenceDiagonalize",

    --WittDecompStuff.m2
    
    --GrothendieckWittClasses.m2    
    "GrothendieckWittClass",
    "baseField",
    "gwClass",
    "gwAdd",
    "gwMultiply",
    
    --BuildingForms.m2
    "diagonalClass",
    "hyperbolicForm",
    "PfisterForm",
    
    --SimplifiedRepresentatives.m2
    "diagonalForm",
    "diagonalEntries",
    "integralDiagonalRep",
    
    --HilbertSymbols.m2
    "HilbertSymbol",
    
    --GWInvariants.m2
    "signature",
    "integralDiscriminant",
    "relevantPrimes",
    "HasseWittInvariant",

    --LocalGlobalDegrees.m2
    "globalA1Degree",
    "localA1Degree",
    
    --IsomorphismOfForms.m2
    "isIsomorphicFormQTommy",
    "isIsomorphicFormQ",
    "gwIsomorphic",
    
    --Isotropy.m2
    "isIsotropicQp",
    "isIsotropic",
    "isAnisotropic",

    --AnisotropicDimension.m2
    "isHyperbolicQp",
    "anisotropicDimensionQp",
    "anisotropicDimensionQQ",
    
    --Decomposition.m2
    "sumDecomposition",
    "sumDecompositionString"
    }

-- Basic arithmetic, p-adic, and commutative algebra operations we will use
load "./A1-Brouwer/ArithmeticMethods.m2"

-- Basic manipulations of matrices we will use
load "./A1-Brouwer/MatrixMethods.m2"

-- For decomposing forms
-- can probably delete this later
load "./A1-Brouwer/WittDecompStuff.m2"

-- Establishing the GrothendieckWittClass type and some basic manipulations
load "./A1-Brouwer/GrothendieckWittClasses.m2"

-- For building new symmetric bilinear forms
load "./A1-Brouwer/BuildingForms.m2"

-- For providing simplified representatives of symmetric bilinear forms
load "./A1-Brouwer/SimplifiedRepresentatives.m2"

-- For Hilbert symbols over p-adic numbers
load "./A1-Brouwer/HilbertSymbols.m2"

-- Invariants of symmetric bilinear forms
load "./A1-Brouwer/GWInvariants.m2"
    
-- Local and global A1-brouwer degrees
load "./A1-Brouwer/LocalGlobalDegrees.m2"

-- Checking if forms are isomorphic
load "./A1-Brouwer/IsomorphismOfForms.m2"

-- For verifying (an)isotropy
load "./A1-Brouwer/Isotropy.m2"

-- Anisotropic dimension
load "./A1-Brouwer/AnisotropicDimension.m2"

-- Finally, decomposing forms
load "./A1-Brouwer/Decomposition.m2"








----------------------------
----------------------------
-- DOCUMENTATION
----------------------------
----------------------------

beginDocumentation()

document{
    Key => A1BrouwerDegrees,
    Headline => "A package for running A1-Brouwer degree computations in Macaulay2",
    }

undocumented {
    }

document {
    Key => {(squarefreePart, QQ), (squarefreePart, ZZ), squarefreePart},
	Headline => "smallest magnitude representative of a square class over the rationals or integers",
	Usage => "squarefreePart(q)",
	Inputs => {
	    QQ => "q" => {"a rational number"},
	    -- ZZ => "n" => {"an integer"}
	    },
	Outputs => {
	    ZZ => { "the smallest magnitude integer in the square class of ", TT "n"}
	    },
	PARA {"Given a rational number (or integer), ", TT "q", ", this command outputs the smallest magnitude integer, ",
                TEX///$m$///, ", such that ", TEX///$q=lm$///, " for some rational number (or integer) ",
				TEX///$l$///, "."},
	EXAMPLE lines ///
		 squarefreePart(15/72)
		 squarefreePart(-1/3)
	 	 ///,
        }

document{
    Key => GrothendieckWittClass,
    Headline => "a new type, intended to capture the isomorphism class of an element of the Grothendieck-Witt ring of a base field",
    PARA {"A ", TT "GrothendieckWittClass" ," object is a type of ", TO2(HashTable, "HashTable"), " encoding the isomorphism class of a non-degenerate symmetric bilinear form ", TEX///$V \times V \to k$///, " over a field ", TEX///$k$///, "."},
    PARA{"Given any basis ", TEX///$e_1,\ldots,e_n$///, " for ", TEX///$V$///, " as a ", TEX///$k$///, "-vector space, we can encode the symmetric bilinear form ", TEX///$\beta$///, " by how it acts on basis elements. That is, we can produce a matrix ", TEX///$\left(\beta(e_i,e_j)\right)_{i,j}$///, ". This is called a ", EM "Gram matrix", " for the symmetric bilinear form. A change of basis will produce a congruent Gram matrix, thus a matrix represents a symmetric bilinear form uniquely up to matrix congruence."},
	
	
    PARA{"A GrothendieckWittClass object can be built from a symmetric ", TO2(matrix, "matrix"), " over a field using the ", TO2(gwClass,"gwClass"), " method."},
    EXAMPLE lines///
    beta = gwClass(matrix(QQ,{{0,1},{1,0}}))
    class beta
    ///,
    PARA{"The underlying matrix representative of a form can be recovered via the ", TT "matrix", " command, and its underlying field can be recovered using ", TO2(baseField,"baseField"), "."},
    EXAMPLE lines///
    beta.matrix
    baseField(beta)
    ///,
    PARA{"For computational purposes, it is often desirable to diagonalize a Gram matrix. Any symmetric bilinear form admits a diagonal Gram matrix representative by ", EM "Sylvester's law of inertia", ", and this is implemented via the ", TO2(diagonalForm, "diagonalForm"), " method."},
    EXAMPLE lines///
    diagonalForm(beta)
    ///,
    PARA{"Once a form has been diagonalized, it is recorded in the cache for ", TT "GrothendieckWittClass", " and can therefore be quickly recovered."},
    EXAMPLE lines///
    beta.cache.diagonalForm
    ///,
    SeeAlso => {"gwClass","diagonalForm","baseField"},
    }

document {
    Key => {(gwClass, Matrix), gwClass},
	Headline => "Grothendieck Witt class of a symmetric matrix",
	Usage => "gwClass(M)",
	Inputs => {
	    Matrix => "M" => {"a symmetric matrix defined over an arbitrary field"}
	    },
	Outputs => {
	    GrothendieckWittClass => { "the isomorphism class of a symmetric bilinear form represented by ", TEX/// $M$///, "." }
	    },
	PARA {"Given a symmetric matrix, ", TEX///$M$///, ", this command outputs an object of type ", TT "GrothendieckWittClass", ". ",
                "This output has the representing matrix, ", TEX///$M$///, ", and the base field of the matrix stored in its CacheTable."},
	EXAMPLE lines ///
		 M := matrix(QQ,{{0,0,1},{0,1,0},{1,0,0}});
		 beta = gwClass(M)
	 	 ///,
	PARA{"The matrix representing a ", TT "GrothendieckWittClass", " element can be recovered using the ", TT "matrix", " command:"},
	EXAMPLE lines ///
	    	beta.matrix
		///,
        PARA{"The base field which the form ", TT "beta", " is implicitly defined over can be recovered with the ", TT "baseField", " method."},
	EXAMPLE lines ///
	    	baseField beta
		///,
		
	SeeAlso => {"baseField","GrothendieckWittClass"}
        }

document {
	Key => {(congruenceDiagonalize, Matrix), congruenceDiagonalize},
	Headline => "diagonalizing a symmetric matrix via congruence",
	Usage => "congruenceDiagonalize(M)",
	Inputs => {
	    Matrix => "M" => {"a symmetric matrix over any field"}
	    },
	Outputs => {
	    Matrix => { "a diagonal matrix congruent to ", TT "M" }
	    },
	PARA {"Given a symmetric matrix ", TEX///$M$///, " over any field, this command gives a diagonal matrix congruent to ", TEX///$M$///,". Note that the order in which the diagonal terms appear is not specified."},
	EXAMPLE lines ///
		 M=matrix(GF(17), {{7, 9}, {9, 6}});
		 congruenceDiagonalize(M)
	 	 ///,
	SeeAlso => {"diagonalForm"}
     	}

document {
    Key => {(baseField, GrothendieckWittClass), baseField},
	Headline => "base field of a Grothendieck Witt class",
	Usage => "baseField(beta)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"the isomorphism class of a symmetric bilinear form"}
	    },
	Outputs => {
	    Ring => { "the base field of the Grothendieck-Witt class ", TT "beta" }
	    },
	PARA {"Given the isomorphism class of a symmetric bilinear form, ", TT "beta", 
                ", this command outputs the base field of the form."},
	EXAMPLE lines ///
		 beta = gwClass(matrix(QQ,{{0,2},{2,0}}));
		 baseField beta
	 	 ///,
    SeeAlso => {"GrothendieckWittClass"}
        }

document {
    Key => {(localAlgebraBasis, List, Ideal), localAlgebraBasis},
	Headline => "produces a basis for a local finitely generated algebra over a field k",
	Usage => "localAlgebraBasis(L,p)",
	Inputs => {
	    List => "L" => {"list of polynomials ", TEX///$f=(f_1, \dots ,f_n)$///, " over the same ring"},
	    Ideal => "p" => {"prime ideal of an isolated zero"}
	    },
	Outputs => {
	    List => {"a list of basis elements of the local k-algebra ", TEX///$Q_p(f)$/// }
	    },
	PARA {"Given an endomorphism of affine space, ", TEX///$f=(f_1,\dots ,f_n)$///,
			", given as a list of polynomials called ", TT "L", " and the prime ideal of an isolated zero, this command returns a list of basis elements of the local k-algebra ", TEX///$Q_p(f)$///, " by computing a normal basis for ", TEX///$(I:(I:p^{\infty}))}$///, " (vis. [S02, Proposition 2.5])."},
	EXAMPLE lines ///
		 QQ[x,y];
		 f = {x^2+1-y,y};
		 p = ideal(x^2+1,y);
		 localAlgebraBasis(f,p) 
	 	 ///,
        PARA{EM "Citations:"},
    UL{
	
	{"[S02] B. Sturmfels, ", EM "Solving Systems of Polynomial Equations,", " American Mathematical Society, 2002."},
	},
    SeeAlso => {"localA1Degree"}
}

document {
    Key => {(sumDecomposition, GrothendieckWittClass), sumDecomposition},
    Headline => "produces a simplified diagonal representative of a Grothendieck Witt class",
    Usage => "sumDecomposition(beta)",
    Inputs => {
        GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over a field ", TEX///$k$///, "."},
    },
    Outputs => { 
	GrothendieckWittClass => {"a diagonal representative of the Grothendieck Witt class of the input form"},
	--String => {"The decomposition as a sum of hyperbolic and rank one forms."},
	},
    PARA {"Given a symmetric bilinear form ", TT"beta", " over a field ", TEX///$k$///, ", we decompose it as a sum of some number of hyperbolic and rank one forms."},
    EXAMPLE lines ///
    M = matrix(RR,{{2.091,2.728,6.747},{2.728,7.329,6.257},{6.747,6.257,0.294}});
    beta = gwClass(M);
    sumDecomposition(beta)
    ///,
    PARA {"Over ", TEX///$\mathbb{R}$///, " there are only two square classes and a form is determined uniquely by its rank and signature [L05, II Proposition 3.2]. A form defined by the ", TEX///$3\times 3$///, " Gram matrix ", TT"M", " above is isomorphic to the form ", TEX///$\langle 1,-1,1\rangle $///, "."},
    EXAMPLE lines ///
    M = matrix(GF(13),{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
    beta = gwClass(M);
    sumDecomposition(beta)
    ///,
    PARA{EM "Citations:"},
    UL{
	
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	},
    SeeAlso => {"sumDecompositionString"},
}

document {
    Key => {(sumDecompositionString, GrothendieckWittClass), sumDecompositionString},
    Headline => "produces a simplified diagonal representative of a Grothendieck Witt class",
    Usage => "sumDecompositionString(beta)",
    Inputs => {
        GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over a field ", TEX///$k$///, "."},
    },
    Outputs => { 
	--GrothendieckWittClass => {"a diagonal representative of the Grothendieck Witt class of the input form"},
	String => {"The decomposition as a sum of hyperbolic and rank one forms."},
	},
    PARA {"Given a symmetric bilinear form ", TT"beta", " over a field ", TEX///$k$///, ", we return a simplified diagonal form of ", TT"beta","."},
    EXAMPLE lines ///
    M = matrix(RR,{{2.091,2.728,6.747},{2.728,7.329,6.257},{6.747,6.257,0.294}});
    beta = gwClass(M);
    sumDecompositionString(beta)
    ///,
    PARA {"Over ", TEX///$\mathbb{R}$///, " there are only two square classes and a form is determined uniquely by its rank and signature [L05, II Proposition 3.2]. A form defined by the ", TEX///$3\times 3$///, " Gram matrix ", TT"M", " above is isomorphic to the form ", TEX///$\langle 1,-1,1\rangle $///, "."},
    EXAMPLE lines ///
    M = matrix(GF(13),{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
    beta = gwClass(M);
    sumDecompositionString(beta)
    ///,
    PARA {"Over ", TEX///$\mathbb{F}_{q}$///, " forms can similarly be diagonalized. In this case as ", TEX///$\langle 1,-1,1,-6 \rangle$///, "."},
    PARA{EM "Citations:"},
    UL{
	
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	},
    SeeAlso => {"sumDecomposition"},
}

document {
    Key => {(globalA1Degree, List), globalA1Degree},
    Headline => "Computes a global A1-Brouwer degree of a list of n polynomials in n variables over a field k",
    Usage => "globalA1Degree(L)",
    Inputs => {
	List => "L" => {"list of polynomials ", TEX///$f = (f_1, \ldots, f_n)$///, " in the polynomial ring ", TEX///$k[x_1,\ldots,x_n]$///, " over a field ", TEX///$k$///, "."}
	},
    Outputs => {
	GrothendieckWittClass => {"The class ", TEX///$\text{deg}^{\mathbb{A}^1}(f)$///, " in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."}
	},
    PARA{"Given an endomorphism of affine space ", TEX///$f=(f_1,\dots ,f_n) \colon \mathbb{A}^n_k \to \mathbb{A}^n_k$///, " with isolated zeros, we may compute its ", TEX///$\mathbb{A}^1$///, EM "-Brouwer degree", " valued in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."
	},
    PARA{"The ",
	TEX///$\mathbb{A}^1$///,
	EM "-Brouwer degree",
	" first defined by Morel [M12] is an algebrao-geometric enrichment of the classical topological Brouwer degree. Using the tools of motivic homotopy theory, one may associate to an endomorphism of affine space the isomorphism class of a symmetric bilinear form whose invariants encode geometric data about how the morphism transforms space."
        },
    PARA{"Such an association appears in the work of Eisenbud-Levine [EL77] and Khimshiashvili [K77], wherein the authors develop a symmetric bilinear form whose signature computes the local degree of a smooth map of real manifolds in the case where the Jacobian may vanish on an affine chart. This was proven to agree with Morel's ", TEX///$\mathbb{A}^1$///, "-Brouwer degree in work of Kass and Wickelgren [KW19]. A similar production of a symmetric bilinear form is given by work of Scheja and Storch [SS76], which develops a symmetric bilinear form attached to a complete intersection. This was also shown to align with the ", TEX///$\mathbb{A}^1$///, "-Brouwer degree in [BW23]."},
    PARA{"Following recent work of B. McKean and Pauli [BMP23], the ", TEX///$\mathbb{A}^1$///, "-Brouwer degree can be computed as a multivariate ", EM "Bezoutian bilinear form.", " The algorithms for producing such a form are developed here."},
    EXAMPLE lines ///
    QQ[x];
    f = {x^2+1};
    globalA1Degree(f)
    ///,
    PARA{"The previous example produces a rank two form with signature zero. This corresponds to the fact that the degree of the complex map ", TEX///$\mathbb{C}\to\mathbb{C},\ z\mapsto z^2$///, " has degree two, while the associated real map ", TEX///$\mathbb{R}\to\mathbb{R},\ x\mapsto x^2$///, " has global degree zero."},
    PARA{"Following [M21] we may think about the ",TEX///$\mathbb{A}^1$///, "-Brouwer degree ", TEX///$\text{deg}^{\mathbb{A}^1}(f)$///, " as a quadratically enriched intersection multiplicity of the hyperplanes ", TEX///$V(f_1)\cap \cdots \cap V(f_n).$///, " As a toy example, consider the curve ", TEX///$y=x(x-1)(x+1)$///, " intersecting the ", TEX///$x$///, "-axis."},
    EXAMPLE lines ///
    QQ[x,y];
    f = {x^3 - x^2 - y, y};
    globalA1Degree(f)
    ///,
    PARA{"The rank of this form is three, as cubics over the complex numbers have three roots counted with multiplicity. This form has signature one, which indicates that when the cubic intersects the ", TEX///$x$///, "-axis, when the three points of intersection are counted with a sign corresponding to a right hand rule, the sum equals one."},
    
    PARA{"The global ", TEX///$\mathbb{A}^1$///, "-Brouwer degree can be computed as a sum over the ", TO2(localA1Degree,"local degrees"), " at the points in the zero locus of the morphism. In the previous example, we see that ", TEX///$V(f)$///, " consists of three points on the affine line. We can compute local degrees at all of these:"},
    EXAMPLE lines///
    QQ[x,y];
    f = {x^3 - x^2 - y, y};
    point1 = ideal{x-1, y};
    localA1Degree(f,point1)
    ///,
    
    PARA{EM "Citations:"},
    UL{
	
	{"[BW23] T. Bachmann, K. Wickelgren, ", EM "Euler classes: six-functors formalism, dualities, integrality and linear subspaces of complete intersections,", " J. Inst. Math. Jussieu, 2023."},
	{"[EL77] D. Eisenbud, H. Levine, ", EM "An algebraic formula for the degree of a C^infinity map germ,", " Annals of Mathematics, 1977."},
	{"[K77] G. Khimshiashvili, ", EM "The local degree of a smooth mapping,", " Sakharth. SSR Mcn. Akad. Moambe, 1977."},
	{"[KW19] J. Kass, K. Wickelgren, ", EM "The class of Eisenbud-Khimshashvili-Levine is the local A1-Brouwer degree,", " Duke Math J., 2019."},
	{"[M21] S. McKean, ", EM "An arithmetic enrichment of Bezout's Theorem,", " Math. Ann., 2021."},
	{"[M12] F. Morel, ", EM "A1-Algebraic topology over a field,", " Springer Lecture Notes in Mathematics, 2012."},
	{"[BMP23] T. Brazelton, S. McKean, S. Pauli, ", EM "Bezoutians and the A1-Degree,", " Algebra & Number Theory, 2023."},
	{"[SS76] S. Scheja, S. Storch, ", EM "Uber Spurfunktionen bei vollstandigen Durchschnitten,", " J. Reine Angew. Math., 1975."},
	},
    SeeAlso => {"localA1Degree", "sumDecomposition", "sumDecompositionString"}
    }

document {
    Key => {(localA1Degree, List, Ideal), localA1Degree},
    Headline => "Computes a local A1-Brouwer degree of a list of n polynomials in n variables over a field k at a prime ideal in the zero locus",
    Usage => "locallA1Degree(L,p)",
    Inputs => {
	List => "L" => {"list of polynomials ", TEX///$f = (f_1, \ldots, f_n)$///, " in the polynomial ring ", TEX///$k[x_1,\ldots,x_n]$///, " over a field ", TEX///$k$///, "."},
	Ideal => "p" => {"a prime ideal ", TEX///$p \trianglelefteq k[x_1,\ldots,x_n]$///, " in the zero locus ", TEX///$V(f)$///, "."},
	},
    Outputs => {
	GrothendieckWittClass => {"The class ", TEX///$\text{deg}_p^{\mathbb{A}^1}(f)$///, " in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."}
	},
    PARA{"Given an endomorphism of affine space ", TEX///$f=(f_1,\dots ,f_n) \colon \mathbb{A}^n_k \to \mathbb{A}^n_k$///, " and an isolated zero ", TEX///$p\in V(f)$/// ,", we may compute its local ", TEX///$\mathbb{A}^1$///, EM "-Brouwer degree", " valued in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///, "."
	},
    PARA{"For historical and mathematical background, see ", TO2(globalA1Degree, "global A1-degrees"), "."},
    EXAMPLE lines ///
    T1 = QQ[z_1..z_2];
    f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
    f1GD = globalA1Degree(f1);
    q=ideal {z_1,z_2};
    r=ideal {z_1-1,z_2^2-(9/85)};
    f1LDq= localA1Degree(f1,q)
    f1LDr= localA1Degree(f1,r)
    f1LDsum = gwAdd(f1LDq, f1LDr)
    ///,
    PARA{"The sum of the local A1-degrees is equal to the global A1-degree:"},
    EXAMPLE lines///
    gwIsomorphic(f1GD,f1LDsum)
    ///,
    SeeAlso => {"globalA1Degree", "sumDecomposition", "sumDecompositionString"}
    }

document{
    Key => {(signature, GrothendieckWittClass), signature},
    Headline => "Outputs the signature of a symmetric bilinear form over the real or rational numbers",
    Usage => "signature(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"A symmetric bilinear form defined over ", TEX///$\mathbb{Q}$///, " or ", TEX///$\mathbb{R}$///, "."},
	},
    Outputs => {
	ZZ => "n" => {"The ", EM "signature", " of the symmetric bilinear form ", TEX///$\beta$///, "."},
	},
    PARA{"Given a symmetric bilinear form, after diagonalizing it, we can consider the number of positive entries minus the number of negative entries appearing along the diagonal. This is the ", EM "signature", " of a symmetric bilinear form, and is one of the primary invariants we use to classify forms. For more information see ", TO2(gwIsomorphic,"gwIsomorphic"), "."},
    EXAMPLE lines ///
    M = matrix(RR,{{0,0,1},{0,1,0},{1,0,0}});
    beta = gwClass(M);
    signature(beta)
    ///,
    SeeAlso => {"gwIsomorphic", "sumDecomposition", "sumDecompositionString"}
    }

document{
    Key => {(gwIsomorphic, GrothendieckWittClass, GrothendieckWittClass), gwIsomorphic},
    Headline => "Determines whether two Grothendieck Witt classes are isomorphic over CC, RR, QQ, or a finite field.",
    Usage => "gwIsomorphic(alpha,beta)",
    Inputs => {
	GrothendieckWittClass => "alpha" => {"Any Grothendieck-Witt class ", TEX///$\alpha$///, "."},
	GrothendieckWittClass => "beta" => {"Any Grothendieck-Witt class ", TEX///$\beta$///, "."},
	},
    Outputs => {
	Boolean => {"returns true or false depending on whether two Grothendieck Witt classes are equal in the Grothendieck-Witt ring"},
	},
    PARA{"Given two matrices representing symmetric bilinear forms over a field ", TEX///$k$///, ", it is a fundamental question to ask when they are representing the same symmetric ibilinear form, i.e. when they are equal in the Grothendieck-Witt ring ", TEX///$\text{GW}(k)$///,"."},
    
    PARA{EM "Sylvester's Law of Inertia", " proves that any symmetric bilinear form can be diagonalized into a block sum of rank one symmetric bilinear forms. Since the rank one forms ", TEX///$\langle a \rangle \colon k \times k \to k$///, ", ", TEX///$(x,y) \mapsto axy$///, " and ", TEX///$\langle ab^2 \rangle \colon k \times k \to k$///, ", ", TEX///$(x,y) \mapsto ab^2xy$///, " differ by a change of basis in the ground field, it follows they are isomorphic (provided that ", TEX///$a,b\ne 0$///, "). Thus after diagonalizing a form, it suffices to consider the square class of each entry appearing along the diagonal. Consider the following example."},
    EXAMPLE lines ///
    alpha = gwClass(matrix(CC,{{2,3,1},{3,-1,0},{1,0,0}}))
    beta = gwClass(matrix(CC,{{2,4,-1},{4,5,7},{-1,7,9}}))
    gwIsomorphic(alpha,beta)
    ///,
    PARA{"The two forms are isomorphic since they can be diagonalized, after which they can be rewritten as the identity matrix after a chance of basis, since every nonzero element is a square class over ", TEX///$\mathbb{C}$///, " (the same is true for any quadratically closed field). Thus we have that the ", EM "rank", " of a form completely determines it over the complex numbers. That is, it provides an isomorphism ", TEX///$\text{GW}(\mathbb{C}) \to \mathbb{Z}$/// ,"."},
    PARA{"Over the reals, the story is a bit different. Since there are two classes over the reals, ", TEX///$ \mathbb{R}^\times / \left(\mathbb{R}^\times\right)^2 \cong \left\{\pm 1\right\}$///," we have a further invariant which classifies symmetric bilinear forms, called the ", TO2(signature, "signature"), ". This is computed as first diagonalizing, then taking the number of positive entries appearing on the diagonal minus the number of negative entries appearing on the diagonal."},
    EXAMPLE lines ///
    gamma = gwClass(matrix(RR,{{1,0,0},{0,-1,0},{0,0,1}}));
    signature(gamma)
    ///,
    PARA{"Rank and signature completely classify symmetric bilinear forms over the reals."},
    EXAMPLE lines ///
    delta = gwClass(matrix(RR,{{0,0,1},{0,1,0},{1,0,0}}));
    gwIsomorphic(gamma,delta)
    ///,
    PARA{"Over finite fields, rank is still an invariant of a form, however signature no longer makes sense as the field is not totally ordered. Instead we consider the ", EM "discriminant", " of the non-degenerate symmetric bilinear form, which is the determinant of any Gram matrix representing the form. The discriminant is well-defined once we consider its target as landing in square classes of the field. ", TEX///$\text{GW}(\mathbb{F}_q) \to \mathbb{F}_q^\times / \left(\mathbb{F}_q^\times\right)^2$///, ". Over finite fields we look to compare both rank and discriminant of a form."},
    EXAMPLE lines///
    alphaF = gwClass(matrix(GF(7),{{1,2,2},{2,0,1},{2,1,5}}))
    betaF = gwClass(matrix(GF(7),{{2,5,1},{5,6,1},{1,1,3}}))
    gammaF = gwClass(matrix(GF(7),{{0,2,4},{2,3,3},{4,3,1}}))
    det(alphaF.matrix)    
    det(betaF.matrix)
    det(gammaF.matrix)
    ///,
    PARA{"We see that ", TEX///$\text{disc}(\alpha)$///, " is a square, while the discriminants of ", TEX///$\beta$///, " and ", TEX///$\gamma$///, " are not. Therefore we see that ", TEX///$\beta \cong \gamma$///, " but neither of them are isomorphic to ", TEX///$\alpha$///, "."},
    EXAMPLE lines///
    gwIsomorphic(alphaF,betaF)
    gwIsomorphic(alphaF,gammaF)
    gwIsomorphic(betaF,gammaF)
    ///,
    PARA{"Over the rationals, further invariants must be considered. We first check if the rank, discriminant, and signature (when considered as a real form) all agree. If so, we must further check whether the ", EM "Hasse-Witt invariants", " agree at all primes. This is an instance of the ", EM "Hasse-Minkowski principle", " which states that quadratic forms are isomorphic over a global field if and they are isomorphic over all its completions (see [S73, IV Theorem 7] or [L05, VI.3.3])."},
    PARA{"The ", EM "Hasse-Witt invariant", " of a diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " over a field ", TEX///$K$///, " is defined to be the product ", TEX///$\prod_{i<j} \left( \phi(a_i,a_j) \right)$///, " where ", TEX///$\phi \colon K \times K \to \left\{\pm 1\right\}$///, " is any ", EM "symbol", " (see e.g. [MH73, III.5.4] for a definition). It is a classical result of Hilbert that over a local field of characteristic not equal to two, there is one and only symbol, ", TEX///$(-,-)_p$///,  " called the ", EM "Hilbert symbol", " ([S73, Chapter III]) computed as follows:"},
    PARA{TEX///$(a,b)_p = \begin{cases} 1 & z^2 = ax^2 + by^2 \text{ has a nonzero solution in } K^3 \\ -1 & \text{otherwise.} \end{cases}$///},
    PARA{"Consider the following example, where we observe that ", TEX///$z^2 = 2x^2 + y^2$///," does admit nonzero solutions mod 7, in particular ", TEX///$(x,y,z) = (1,0,3)$///, ":"},
    EXAMPLE lines///
    HilbertSymbol(2,1,7)
    ///,
    PARA{"The Hasse invariant will be 1 for almost all primes. In particular after diagonalizing a form ", TEX///$\beta \cong \left\langle a_1,\ldots,a_n\right\rangle$///, " then the Hasse invariant at a prime ", TEX///$p$///, " will be 1 automatically if ", TEX///$p\nmid a_i$///, " for all ", TEX///$i$///, ". Thus we only have finitely many Hasse invariants to compare for any pair of symmetric bilinear forms."},
    EXAMPLE lines///
    alphaQ = gwClass(matrix(QQ,{{1,4,7},{4,3,2},{7,2,-1}}))
    betaQ = gwClass(matrix(QQ,{{0,0,1},{0,2,7},{1,7,3}}))
    gwIsomorphic(alphaQ,betaQ)
    ///,
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	{"[MH73] Milnor and Husemoller, ", EM "Symmetric bilinear forms,", " Springer-Verlag, 1973."},
    },
    SeeAlso => {"signature", "sumDecomposition", "sumDecompositionString"}
}

document{
    Key => {(HilbertSymbol, ZZ,ZZ,ZZ), HilbertSymbol},
    Headline => "Computes the Hilbert symbol of two integers at a prime",
    Usage => "HilbertSymbol(a,b,p)",
    Inputs => {
	ZZ => "a" => {"Any integer, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	ZZ => "b" => {"Any integer, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	ZZ => "p" => {"Any prime number."},
	},
    Outputs => {
	ZZ => {"The ", EM "Hilbert symbol ", TEX///$(a,b)_p$///, "."},
	},
    PARA{"The ", EM "Hasse-Witt invariant", " of a diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " over a field ", TEX///$K$///, " is defined to be the product ", TEX///$\prod_{i<j}  \phi(a_i,a_j)$///, " where ", TEX///$\phi \colon K \times K \to \left\{\pm 1\right\}$///, " is any ", EM "symbol", " (see e.g. [MH73, III.5.4] for a definition). It is a classical result of Hilbert that over a local field of characteristic not equal to two, there is one and only symbol, ", TEX///$(-,-)_p$///,  " called the ", EM "Hilbert symbol", " ([S73, Chapter III]) computed as follows:"},
    PARA{TEX///$(a,b)_p = \begin{cases} 1 & z^2 = ax^2 + by^2 \text{ has a nonzero solution in } K^3 \\ -1 & \text{otherwise.} \end{cases}$///},
    PARA{"Consider the following example, where we observe that ", TEX///$z^2 = 2x^2 + y^2$///," does admit nonzero solutions mod 7, in particular ", TEX///$(x,y,z) = (1,0,3)$///, ":"},
    EXAMPLE lines///
    HilbertSymbol(2,1,7)
    ///,
    PARA{"Computing Hasse-Witt invariants is a key step in classifying symmetric bilinear forms over the rational numbers, and in particular certifying their ", TO2(isIsotropic, "(an)isotropy"), "."},
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[MH73] Milnor and Husemoller, ", EM "Symmetric bilinear forms,", " Springer-Verlag, 1973."},
    },
}

document {
    Key => {(diagonalForm, GrothendieckWittClass), diagonalForm},
	Headline => "produces a diagonalized form for any Grothendieck-Witt class",
	Usage => "diagonalForm(beta)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"any class in ", TEX///$\text{GW}(k)$///," where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."}
	    },
	Outputs => {
	    GrothendieckWittClass => {"a form isomorphic to ", TEX///$\beta$///, " with a diagonal Gram matrix"}
	    },
	PARA {"Given a symmetric bilinear form, this method calls the ", TO2(congruenceDiagonalize,"congruenceDiagonalize"), " command in order to produce a diagonal symmetric bilinear form isomorphic to ", TEX///$\beta$///, "."},
	EXAMPLE lines ///
	beta = gwClass(matrix(QQ,{{0,0,2},{0,2,0},{2,0,0}}));
	diagonalForm(beta)
        ///,
	PARA{"Note that the  ", TO2(GrothendieckWittClass, "GrothendieckWittClass"), " type caches diagonal versions of a form once they've been computed. We can recover this quickly in the following way."},
	EXAMPLE lines///
	beta.cache.diagonalForm
	///,
	SeeAlso => {"congruenceDiagonalize","diagonalEntries"}
	}

document {
    Key => {(diagonalEntries, GrothendieckWittClass), diagonalEntries},
	Headline => "extracts a list of diagonal entries for a GrothendieckWittClass",
	Usage => "diagonalEntries(beta)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"any class in ", TEX///$\text{GW}(k)$///," where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."}
	    },
	Outputs => {
	    List => "L" => {"a list ", TEX///$L = \{a_1,\ldots,a_n\}$///, " of elements ", TEX///$a_i\in k$///, " so that ", TEX///$\beta \cong \langle a_1,\ldots,a_n\rangle$///, "."}
	    },
	PARA {"Given a diagonal form, ", TT "diagonalEntries", " will extract the elements along the diagonal."},
	EXAMPLE lines ///
	beta = gwClass(matrix(QQ,{{3,0,0},{0,2,0},{0,0,7}}))
	diagonalEntries beta
        ///,
	PARA{"If the form is not given with a diagonal representative, this method will first diagonalize it."},
	EXAMPLE lines///
	gamma = gwClass(matrix(RR,{{0,0,1},{0,1,0},{1,0,0}}))
	diagonalEntries gamma
	///,
	SeeAlso => {"diagonalForm", "congruenceDiagonalize"}
	}
    
    
document {
    Key => {(integralDiagonalRep, GrothendieckWittClass), integralDiagonalRep},
	Headline => "given a rational symmetric bilinear form, outputs a diagonal representative with integral entries",
	Usage => "integralDiagonalRep(beta)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"any class in ", TEX///$\text{GW}(\mathbb{Q})$///,"."}
	    },
	Outputs => {
	    GrothendieckWittClass => {"a form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " isomorphic to ", TEX///$\beta$///, " with each ", TEX///$a_i \in \mathbb{Z}$///, "."}
	    },
	EXAMPLE lines ///
	M = matrix(QQ,{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
	beta = gwClass(M);
	integralDiagonalRep(beta)
        ///,
	SeeAlso => {"diagonalForm", "diagonalEntries", "congruenceDiagonalize"}
	}


document{
    Key => {(isIsotropic, GrothendieckWittClass), isIsotropic},
    Headline => "Determines whether a Grothendieck-Witt class is isotropic",
    Usage => "isIsotropic(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(k)$///, " where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."},
	},
    Outputs => {
        Boolean => {"Whether ", TEX///$\beta$///, " is isotropic"},
	},
    PARA{"Recall a symmetric bilinear form ", TEX///$\beta$///, " is said to be ", EM "isotropic", " if there exists a nonzero vector ", TEX///$v$///, " for which ", TEX///$\beta(v,v) = 0$///, ". Witt's decomposition theorem implies that a non-degenerate symmetric bilinear form decomposes uniquely into an isotropic and an anisotropic part. Certifying (an)isotropy is then an important computational problem when working with the Grothendieck-Witt ring."},
    PARA{"Over ", TEX///$\mathbb{C}$///, ", any form of rank two or higher contains a copy of the hyperbolic form, and hence is isotropic. Thus we can determine isotropy simply by a consideration of rank."},
    EXAMPLE lines///
    isIsotropic(gwClass(matrix(CC,{{3}})))
    isIsotropic(gwClass(matrix(CC,{{2,0},{0,5}})))
    ///,
    PARA{"Forms over ", TEX///$\mathbb{R}$///, " are anisotropic if and only if all its diagonal entries are positive or are negative."},
    EXAMPLE lines///
    isIsotropic(gwClass(matrix(RR,{{3,0,0},{0,5,0},{0,0,7}})))
    isIsotropic(gwClass(matrix(RR,{{0,2},{2,0}})))
    ///,
    PARA{"Over finite fields, a form is anisotropic so long as it is nondegenerate, of rank ", TEX///$\le 2$///," and not isomorphic to the hyperbolic form."},
    EXAMPLE lines///
    isIsotropic(gwClass(matrix(GF(7),{{1,0,0},{0,1,0},{0,0,1}})))
    isIsotropic(gwClass(matrix(GF(7),{{3,0},{0,3}})))
    ///,
    PARA{"Over ", TEX///$\mathbb{Q}$///, " things become a bit more complicated. We can exploit the local-to-global principle for isotropy (the ", EM "Hasse-Minkowski principle", "), which states that a form is isotropic over ", TEX///$\mathbb{Q}$///, " if and only if it is isotropic over all its completions, meaning all the ", TEX///$p$///, "-adic numbers and ", TEX///$\mathbb{R}$///, " [L05, VI.3.1]. We note, however, the classical result that all forms of rank ", TEX///$\ge 5$///, " in ", TEX///$\mathbb{Q}_p$///, " are isotropic [S73, IV Theorem 6]. Thus isotropy in this range of ranks is equivalent to checking it over the real numbers."},
    EXAMPLE lines///
    beta = gwClass(matrix(QQ,{{1, 0, 2, 0, 3}, {0, 6, 1, 1, -1},{2, 1, 5, 2, 0}, {0, 1, 2, 4, -1}, {3, -1, 0,-1, 1}}));
    isIsotropic(beta)
    diagonalForm(beta)
    ///,
    PARA{"For forms of rank ", TEX///$\le 4$///, " we should understand when the form is isotropic over local fields."},
    PARA{"Ternary forms are isotropic away from primes dividing the coefficients of the form in a diagonal basis by e.g. [L05, VI.2.5(2)], so there are only finitely many things to check. Over these relevant primes, isotropy of a form ", TEX///$\beta \in \text{GW}(\mathbb{Q})$///, " over ", TEX///$\mathbb{Q}_p$///," is equivalent to the statement that ", TEX///$(-1,-\text{disc}(\beta))_p = H(\beta)$///, " where ", TEX///$H(\beta)$///, " denotes the Hasse-Witt invariant attached to ", TEX///$\beta$///, " and ", TEX///$(-,-)_p$///," is the ", TO2(HilbertSymbol, "Hilbert Symbol"), "."},
    PARA{"A binary form ", TEX///$q$///, " is isotropic if and only if it is isomorphic to the hyperbolic form, which implies in particular that the rank, signature, and discriminant of ", TEX///$q$///, " agree with that of ", TEX///$\mathbb{H}=\langle 1,-1\rangle$///, ". " },
    PARA{"TODO --- rewrite this whole documentation using anisotropicDimension to compute isIsotropic / isAnisotropic booleans"},
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
    },
    SeeAlso => {"isAnisotropic", "HilbertSymbol", "signature"}
}

document{
    Key => {(isAnisotropic, GrothendieckWittClass), isAnisotropic},
    Headline => "Determines whether a Grothendieck-Witt class is anisotropic",
    Usage => "isAnisotropic(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(k)$///, " where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."},
	},
    Outputs => {
        Boolean => {"Whether ", TEX///$\beta$///, " is anisotropic"},
	},
    PARA{"This is the negation of the boolean-valued ", TO2(isIsotropic,"isIsotropic"), ". See documentation there."},
    SeeAlso => {"isIsotropic"}   
}




document{
    Key => {(isIsotropicQp, GrothendieckWittClass,ZZ), isIsotropicQp},
    Headline => "determines whether a rational form is isotropic locally at a prime",
    Usage => "isIsotropicQp(beta,p)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///, "."},
	ZZ => "p" => {"a prime"},
	},
    Outputs => {
        Boolean => {"Whether ", TEX///$\beta$///, " is isotropic over ", TEX///$\mathbb{Q}_p$///,"."},
	},
    PARA{"Every rank one nondegenerate form is anisotropic, while every form of rank ", TEX///$\ge 5$///, " is isotropic."},
    EXAMPLE lines ///
    isIsotropicQp(gwClass(matrix(QQ,{{2}})),7)
    ///,
    PARA{"For ranks two, three, and four, we can use various criteria for isotropy as in [S73, IV Theorem 6]."},
    PARA{EM "Citations:"},
    UL{	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},	
      },
    SeeAlso => {"isIsotropic", "isAnisotropic","HilbertSymbol"},
    
}



document {
    Key => {(gwAdd, GrothendieckWittClass, GrothendieckWittClass), gwAdd},
    Headline => "The direct sum of two Grothendieck-Witt classes",     
    Usage => "gwAdd(beta, gamma)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"the isomorphism class of a non-degenerate symmetric bilinear form represented by a matrix ", TT "M"},
	    GrothendieckWittClass => "gamma" => {"the isomorphism class of a non-degenerate symmetric bilinear form represented by a matrix ", TT "N"},
	    }, 
	Outputs => { 
	    GrothendieckWittClass => {"the isomorphism class of the direct sum of the bilinear forms represented by the matrices ", TT "M", " and ", TT "N"},
	    },
	PARA {"This computes the direct sum of the Grothendieck-Witt classes ",TT "beta"," and ",TT "gamma","."},
	EXAMPLE lines ///
		 M = matrix(QQ,{{1,0},{0,1}});
		 N = matrix(QQ, {{1, 2}, {2, 5}});
		 beta = gwClass(M);
		 gamma = gwClass(N);
    	    	 gwAdd(beta, gamma)
	 	 ///,
    SeeAlso => {"GrothendieckWittClass", "gwClass", "gwMultiply"}
}

document {
    Key => {(gwMultiply, GrothendieckWittClass, GrothendieckWittClass), gwMultiply},
    Headline => "The tensor product of two Grothendieck-Witt classes",     
	Usage => "gwMultiply(beta, gamma)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"the isomorphism class of a non-degenerate symmetric bilinear form represented by a matrix ", TT "M"},
	    GrothendieckWittClass => "gamma" => {"the isomorphism class of a non-degenerate symmetric bilinear form represented by a matrix ", TT "N"},
	    }, 
	Outputs => { 
	    GrothendieckWittClass => {"the isomorphism class of the tensor product of the bilinear forms represented by the matrices ", TT "M", " and ", TT "N"},
	    },
	PARA {"This computes the tensor product of the Grothendieck-Witt classes ",TT "beta"," and ",TT "gamma","."},
	EXAMPLE lines ///
    	    	 M = matrix(QQ,{{1,0},{0,1}});
		 N = matrix(QQ, {{1, 2}, {2, 5}});
		 beta = gwClass(M);
		 gamma = gwClass(N);
    	    	 gwMultiply(beta, gamma)
	 	 ///,
    SeeAlso => {"GrothendieckWittClass", "gwClass", "gwAdd"}
}


document {
    Key => {diagonalClass, (diagonalClass, Ring, RingElement), (diagonalClass,Ring,ZZ),(diagonalClass,Ring,QQ),(diagonalClass,Ring,Sequence), (diagonalClass,InexactFieldFamily,RingElement), (diagonalClass,InexactFieldFamily,ZZ), (diagonalClass,InexactFieldFamily,QQ),(diagonalClass,InexactFieldFamily,Sequence)},
    Headline => "the Grothendieck-Witt class of a diagonal form",     
	Usage => "diagonalClass(k,a)
	          diagonalClass(k,L)",
	Inputs => {
	    Ring => "k" => {"a field"},
	    RingElement => "a" => {"any element ", TEX///$a\in k$///},
	    Sequence => "L" => {"a list of elements ", TEX///$L = (a_1,\ldots,a_n)$///, " with ", TEX///$a_i \in k$///},
	    }, 
	Outputs => { 
	    GrothendieckWittClass => {"the diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle \in \text{GW}(k)$///},
	    },
	PARA {"Given a sequence of elements ", TEX///$a_1,\ldots,a_n \in k$///, " we can form the diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " defined to be the block sum of each of the rank one forms ", TEX///$\langle a_i \rangle \colon k \times k \to k$///, " ", TEX///$(x,y) \mapsto a_i xy$///, "."},
	EXAMPLE lines ///
    	    	 diagonalClass(QQ,(3,5,7))
	 	 ///,
	PARA{"Inputting a ring element, an integer, or a rational instead of a sequence will produce a rank one form instead. For instance:"},
	EXAMPLE lines ///
	diagonalClass(GF(29),5/13)
	diagonalClass(RR,2)
	///,
    SeeAlso => {"diagonalForm", "congruenceDiagonalize", "diagonalEntries"}
}

document {
    Key => {PfisterForm, (PfisterForm, Ring, RingElement), (PfisterForm,Ring,ZZ),(PfisterForm,Ring,QQ),(PfisterForm,Ring,Sequence), (PfisterForm,InexactFieldFamily,RingElement), (PfisterForm,InexactFieldFamily,ZZ), (PfisterForm,InexactFieldFamily,QQ),(PfisterForm,InexactFieldFamily,Sequence)},
    Headline => "the Grothendieck-Witt class of a Pfister form",     
	Usage => "PfisterForm(k,a)
	          PfisterForm(k,L)",
	Inputs => {
	    Ring => "k" => {"a field"},
	    RingElement => "a" => {"any element ", TEX///$a\in k$///},
	    Sequence => "L" => {"a list of elements ", TEX///$L = (a_1,\ldots,a_n)$///, " with ", TEX///$a_i \in k$///},
	    }, 
	Outputs => { 
	    GrothendieckWittClass => {"the Pfister form ", TEX///$\langle\langle a_1,\ldots,a_n\rangle\rangle \in \text{GW}(k)$///},
	    },
	PARA {"Given a sequence of elements ", TEX///$a_1,\ldots,a_n \in k$///, " we can form the Pfister form ", TEX///$\langle\langle a_1,\ldots,a_n\rangle\rangle$///, " defined to be the rank ", TEX///$2^n$///, " form defined as the product ", TEX///$\langle 1, -a_1\rangle \otimes \cdots \otimes \langle 1, -a_n \rangle$///, "."},
	EXAMPLE lines ///
    	    	 PfisterForm(QQ,(2,6))
	 	 ///,
	PARA{"Inputting a ring element, an integer, or a rational instead of a sequence will produce a one-fold Pfister form instead. For instance:"},
	EXAMPLE lines ///
	PfisterForm(GF(13),-2/3)
	PfisterForm(CC,3)
	///
}


document {
    Key => {hyperbolicForm, (hyperbolicForm, Ring), (hyperbolicForm, Ring, ZZ), (hyperbolicForm, InexactFieldFamily), (hyperbolicForm, InexactFieldFamily, ZZ)},
    Headline => "the Grothendieck-Witt class of a hyperbolic form",     
	Usage => "hyperbolicForm(k)
	          hyperbolicForm(k,n)",
	Inputs => {
	    Ring => "k" => {"a field"},
	    ZZ => "n" => {"an even number, giving an optional rank ", TEX///$n$///, " for a totally hyperbolic form"},
	    }, 
	Outputs => { 
	    GrothendieckWittClass => {"the hyperbolic form ", TEX///$\mathbb{H} = \langle 1, -1\rangle \in \text{GW}(k)$///, " or the totally hyperbolic form ", TEX///$\frac{n}{2}\mathbb{H}$///, " if an optional rank is specified"},
	    },
	PARA {"By default outputs the rank two hyperbolic form over the input field ", TEX///$k$///, "."},
	EXAMPLE lines ///
    	    	 hyperbolicForm(GF(7))
	 	 ///,
	PARA{"Specifying a rank yields a copy of sums of the rank two hyperbolic form. Only even rank inputs are accepted."},
	EXAMPLE lines ///
        hyperbolicForm(RR,4)
	///,
    SeeAlso => {"isAnisotropic", "sumDecomposition", "sumDecompositionString"}
}

document{
    Key => {(legendreBoolean, RingElement), legendreBoolean},
    Headline => "Basic Legendre symbol over a finite field",
    Usage => "legendreBoolean(a)",
    Inputs => {
	RingElement => "a" => {"Any element in a finite field ", TEX///$a\in \mathbb{F}_q$///, "."},
	},
    Outputs =>{
	Boolean => {"Whether ", TEX///$a$///, " is a square in ", TEX///$\mathbb{F}_q$///, "."},
	},
    PARA{"Given an element of a finite field, will return a Boolean checking if it is a square."},
    EXAMPLE lines///
    a = sub(-1,GF(5));
    legendreBoolean(a)
    b = sub(-1,GF(7));
    legendreBoolean(b)
    ///,
    }

document{
    Key => {(integralDiscriminant, GrothendieckWittClass), integralDiscriminant},
    Headline => "outputs an integral discriminant for a rational symmetric bilinear form",
    Usage => "integralDiscriminant(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///, "."},
	},
    Outputs => {
        ZZ => {"An integral square class representative of ", TEX///$\text{disc}(\beta)$///, "."},
	},
    EXAMPLE lines ///
    beta = gwClass(matrix(QQ,{{1,4,7},{4,3,-1},{7,-1,5}}));
    integralDiscriminant(beta)
    integralDiagonalRep(beta)
    ///,
}


----------------------------
----------------------------
-- Testing
----------------------------
----------------------------

-- For debugging: remember Macaulay2 starts counting the first test as Test 0

       
------------------
-- TESTING
------------------

-- Diagonal form testing
TEST ///
print("diagonal form testing");
M1=matrix(RR, {{0, 1}, {1, 0}});
G1=gwClass(M1);
M2=diagonalForm(G1);
assert(M2.matrix===matrix(RR, {{1, 0}, {0, -1}}));
///

TEST ///
M3=matrix(CC, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
assert(M4.matrix===matrix(CC, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
-- Ensure the cache is populated
-- assert(G2.cache.?diagonalForm)
///


TEST ///
M3=matrix(QQ, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
assert(M4.matrix===matrix(QQ,{{1, 0, 0}, {0, -2, 0}, {0, 0, 1/2}}));
///

-- gwTypeTest.m2
TEST ///
M = matrix(QQ,{{1,0},{0,1}});
N = matrix(QQ, {{1, 2}, {2, 5}})
beta = gwClass(M);
gamma = gwClass(N);
assert(baseField(beta) === QQ)
assert(beta.matrix === M)
--Operations within GW-classes
A = gwAdd(beta, gamma);
B = gwMultiply(beta, gamma);
assert(A.matrix === matrix(QQ, {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 2}, {0, 0, 2, 5}}));
assert(B.matrix === matrix(QQ, {{1, 2, 0, 0}, {2, 5, 0, 0}, {0, 0, 1, 2}, {0, 0, 2, 5}}));
---non well-defined GW-classes
M'=matrix(ZZ, {{1, 0}, {0, 1}});
N'=matrix(QQ, {{1, 1}, {1, 1}});
theta=gwClass(M');
sigma=gwClass(N');
assert(isWellDefined(theta) === false);
assert(isWellDefined(sigma) === false);
///


-- degreeTesting.m2
TEST ///
T1 = QQ[x]
f = {x^2}
beta = globalA1Degree(f)
gamma = gwClass(matrix(QQ,{{0,1},{1,0}}))
assert(gwIsomorphic(beta,gamma))
///


TEST ///
T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
f1GD = globalA1Degree(f1);
f1GDmat = f1GD.matrix;
assert((WittDecomp(f1GDmat))==(3,matrix(QQ,{{}})));
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
f1LDq= localA1Degree(f1,q);
f1LDr= localA1Degree(f1,r);
f1LDsum = gwAdd(f1LDq, f1LDr);
assert(gwIsomorphic(f1LDsum, f1GD));



T2 = QQ[w];
f2 = {w^4 + w^3 - w^2 - w};
f2GD= globalA1Degree(f2);
f2GDmat = f2GD.matrix;
assert(WittDecomp(f2GDmat)==(2,matrix(QQ,{{}})));

p=ideal {w+1};
f2LDp = localA1Degree(f2, p);
f2LDpmat = f2LDp.matrix;
assert(WittDecomp(f2LDpmat)==(1,matrix(QQ,{{}})));
s=ideal{w-1};
f2LDs = localA1Degree(f2, s);
t=ideal{w};
f2LDt = localA1Degree(f2, t);
f2LDsum = gwAdd(gwAdd(f2LDp, f2LDs),f2LDt);
assert(gwIsomorphic(f2LDsum, f2GD));
///
    


end
