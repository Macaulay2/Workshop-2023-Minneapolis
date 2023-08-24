document{
    Key => {(hilbertSymbol, QQ, QQ,ZZ), hilbertSymbol},
    Headline => "Computes the Hilbert symbol of two integers or rational numbers at a prime",
    Usage => "hilbertSymbol(a,b,p)",
    Inputs => {
	QQ => "a" => {"Any integer or rational number, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	QQ => "b" => {"Any integer or rational number, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	ZZ => "p" => {"Any integer prime number."},
	},
    Outputs => {
	ZZ => {"The ", EM "Hilbert symbol ", TEX///$(a,b)_p$///, "."},
	},
    PARA{"The ", EM "Hasse-Witt invariant", " of a diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " over a field ", TEX///$K$///, " is defined to be the product ", TEX///$\prod_{i<j}  \phi(a_i,a_j)$///, " where ", TEX///$\phi \colon K \times K \to \left\{\pm 1\right\}$///, " is any ", EM "symbol", " (see e.g. [MH73, III.5.4] for a definition). It is a classical result of Hilbert that over a local field of characteristic not equal to two, there is one and only symbol, ", TEX///$(-,-)_p$///,  " called the ", EM "Hilbert symbol", " ([S73, Chapter III]) computed as follows:"},
    PARA{TEX///$(a,b)_p = \begin{cases} 1 & z^2 = ax^2 + by^2 \text{ has a nonzero solution in } K^3 \\ -1 & \text{otherwise.} \end{cases}$///},
    PARA{"Consider the following example, where we observe that ", TEX///$z^2 = 2x^2 + y^2$///," does admit nonzero solutions mod 7, in particular ", TEX///$(x,y,z) = (1,0,3)$///, ":"},
    EXAMPLE lines///
    hilbertSymbol(2,1,7)
    ///,
    PARA{"Computing Hasse-Witt invariants is a key step in classifying symmetric bilinear forms over the rational numbers, and in particular certifying their ", TO2(isIsotropic, "(an)isotropy"), "."},
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[MH73] Milnor and Husemoller, ", EM "Symmetric bilinear forms,", " Springer-Verlag, 1973."},
    },
}

document{
    Key => {(hilbertSymbolReal, QQ, QQ), hilbertSymbolReal},
    Headline => "Computes the Hilbert symbol of two rational numbers over the reals",
    Usage => "hilbertSymbolReal(a,b,p)",
    Inputs => {
	QQ => "a" => {"Any non-zero integer or rational number, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	QQ => "b" => {"Any non-zero integer or rational number, considered as an element of ", TEX///$\mathbb{Q}_p$///, "."},
	},
    Outputs => {
	ZZ => {"The ", EM "Hilbert symbol ", TEX///$(a,b)_{\mathbb{R}}$///, "."},
	},
    PARA{"The ", EM "Hasse-Witt invariant", " of a diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " over a field ", TEX///$K$///, " is defined to be the product ", TEX///$\prod_{i<j}  \phi(a_i,a_j)$///, " where ", TEX///$\phi \colon K \times K \to \left\{\pm 1\right\}$///, " is any ", EM "symbol", " (see e.g. [MH73, III.5.4] for a definition). It is a classical result of Hilbert that over a local field of characteristic not equal to two, there is one and only symbol, ", TEX///$(-,-)_p$///,  " called the ", EM "Hilbert symbol", " ([S73, Chapter III]) computed as follows:"},
    PARA{TEX///$(a,b)_{\mathbb{R}} = \begin{cases} 1 & z^2 = ax^2 + by^2 \text{ has a nonzero solution in } {\mathbb{R}}^3 \\ -1 & \text{otherwise.} \end{cases}$///},
    PARA{TEX///$(a,b)_{\mathbb{R}}$///," will equal 1 unless both ",TEX///$a,\,b$///," are negative."},
    PARA{"Consider the example, that ",TEX///$z^2=-3x^2-2y^2/3$///," does not admit a non-zero solution. Thus:"}, 
     EXAMPLE lines///
    hilbertSymbolReal(-3,-2/3)=-1
    ///,
    PARA{"Computing Hasse-Witt invariants is a key step in classifying symmetric bilinear forms over the rational numbers, and in particular certifying their ", TO2(isIsotropic, "(an)isotropy"), "."},
    PARA{EM "Citations:"},
    UL{
	
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[MH73] Milnor and Husemoller, ", EM "Symmetric bilinear forms,", " Springer-Verlag, 1973."},
    },
}
