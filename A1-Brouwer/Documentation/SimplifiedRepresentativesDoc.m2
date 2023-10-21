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
	SeeAlso => {"congruenceDiagonalize","diagonalEntries","diagonalFormSimplify"}
	}
    

document {
    Key => {(diagonalFormSimplify, GrothendieckWittClass), diagonalFormSimplify},
	Headline => "produces a diagonalized form for any Grothendieck-Witt class, with simplified terms on the diagonal",
	Usage => "diagonalFormSimplify(beta)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"any class in ", TEX///$\text{GW}(k)$///," where ", TEX///$k$///, " is the rationals, reals, complex numbers, or a finite field."}
	    },
	Outputs => {
	    GrothendieckWittClass => {"a form isomorphic to ", TEX///$\beta$///, " with a diagonal Gram matrix"}
	    },
	PARA {"This is the same method as ", TO2(diagonalForm,"diagonalForm")," but will go a step further to reduce square classes appearing on the diagonal entries of a matrix."},
	EXAMPLE lines ///
	beta = gwClass(matrix(QQ,{{0,4},{4,9}}));
	diagonalForm(beta)
	diagonalFormSimplify(beta)
	///,
	SeeAlso => {"diagonalForm"}
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

