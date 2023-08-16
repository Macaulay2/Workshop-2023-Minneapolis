

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
