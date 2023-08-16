document{
    Key => {(anisotropicDimensionQp, GrothendieckWittClass, ZZ), anisotropicDimensionQp},
    Headline => "Returns the anisotropic dimension of a rational symmetric bilinear form over the p-adics",
    Usage => "anisotropicDimensionQp(beta, p)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///, ". "},
	ZZ => "p" => {"A prime number"},
	},
    Outputs => {
        ZZ => {"The rank of the anisotropic part of ", TEX///$\beta$///, " over ", TEX///$\mathbb{Q}_p$///, "."},
	},
    PARA{"This is an implementation of [KC18, Algorithm 8] in Macaulay2, which certifies anisotropy of rational forms over the ", TEX///$p$///,"-adics. Note that any form of rank ", TEX///$\ge 5$///, " is always isotropic, so this method will return 0, 1, 2, 3, or 4."},
    PARA{EM "Citations:"},
    UL{
	
	{"[KC18] P. Koprowski, A. Czogala, ", EM "Computing with quadratic forms over number fields,", " Journal of Symbolic Computation, 2018."},
    },
    SeeAlso => {"anisotropicDimensionQQ", "anisotropicDimension"}   
}

document{
    Key => {(anisotropicDimensionQQ, GrothendieckWittClass), anisotropicDimensionQQ},
    Headline => "Returns the anisotropic dimension of a symmetric bilinear form over the rationals",
    Usage => "anisotropicDimensionQQ(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///, ". "},
	},
    Outputs => {
        ZZ => {"The rank of the anisotropic part of ", TEX///$\beta$///, "."},
	},
    PARA{"This is an implementation of [KC18, Algorithm 9] in Macaulay2. Using ", TO2(anisotropicDimensionQp,"anisotropicDimensionQp"), " and  ", TO2(signature,"signature"), " we can understand the anisotropic dimension of ", TEX///$\beta$///, " at all its relevant completions. The anisotropic dimension over ", TEX///$\mathbb{Q}$///, " is then computed as the maximum of the anisotropic dimensions over all completions."},
    PARA{EM "Citations:"},
    UL{
	
	{"[KC18] P. Koprowski, A. Czogala, ", EM "Computing with quadratic forms over number fields,", " Journal of Symbolic Computation, 2018."},
    },
    SeeAlso => {"anisotropicDimensionQp", "anisotropicDimension"}   
}



document{
    Key => {(anisotropicDimension, GrothendieckWittClass), anisotropicDimension},
    Headline => "Returns the anisotropic dimension of a symmetric bilinear form",
    Usage => "anisotropicDimension(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(k)$///, " where ", TEX///$k$///, " is the complex numbers, reals, rationals, or a finite field."},
	},
    Outputs => {
        ZZ => {"The rank of the anisotropic part of ", TEX///$\beta$///, "."},
	},
    PARA{"By Witt decomposition, any form decomposes uniquely as ", TEX///$\beta \cong k \mathbb{H} \oplus \beta_a$///," where the form ", TEX///$\beta_a$///," is anisotropic. The rank of ", TEX///$\beta_a$///, " is called the ", EM "anisotropic dimension", " of ", TEX///$\beta$///, "."},
    SeeAlso => {"WittIndex", "anisotropicDimensionQp", "anisotropicDimension"}   
}

document{
    Key => {(WittIndex, GrothendieckWittClass), WittIndex},
    Headline => "Returns the Witt index of a symmetric bilinear form",
    Usage => "WittIndex(beta)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"Any class ", TEX///$\beta\in\text{GW}(k)$///, " where ", TEX///$k$///, " is the complex numbers, reals, rationals, or a finite field."},
	},
    Outputs => {
        ZZ => {"The rank of the totally isotropic part of ", TEX///$\beta$///, "."},
	},
    PARA{"By Witt decomposition, any form decomposes uniquely as ", TEX///$\beta \cong k \mathbb{H} \oplus \beta_a$///," where the form ", TEX///$\beta_a$///," is anisotropic. The integer ", TEX///$2k$///, " is called the ", EM "Witt index", " of ", TEX///$\beta$///, "."},
    SeeAlso => {"anisotropicDimension", "anisotropicDimensionQp", "anisotropicDimension"}   
}


