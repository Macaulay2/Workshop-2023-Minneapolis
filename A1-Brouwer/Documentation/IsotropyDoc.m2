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
    PARA{"Ternary forms are isotropic away from primes dividing the coefficients of the form in a diagonal basis by e.g. [L05, VI.2.5(2)], so there are only finitely many things to check. Over these relevant primes, isotropy of a form ", TEX///$\beta \in \text{GW}(\mathbb{Q})$///, " over ", TEX///$\mathbb{Q}_p$///," is equivalent to the statement that ", TEX///$(-1,-\text{disc}(\beta))_p = H(\beta)$///, " where ", TEX///$H(\beta)$///, " denotes the Hasse-Witt invariant attached to ", TEX///$\beta$///, " and ", TEX///$(-,-)_p$///," is the ", TO2(hilbertSymbol, "Hilbert Symbol"), "."},
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
