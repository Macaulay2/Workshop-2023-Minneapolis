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
