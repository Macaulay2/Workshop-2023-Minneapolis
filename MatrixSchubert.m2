newPackage(
    "MatrixSchubert",
    AuxiliaryFiles => true,
    Version => "0.1", 
    Date => "",
    Authors => {
        {Name => "Ayah Almousa", 
            Email => "almou007@umn.edu", 
            HomePage => "http://sites.google.com/view/ayah-almousa"},
	{Name=> "Sean Grate",
	    Email => "sean.grate@auburn.edu",
	    HomePage => "https://seangrate.com/"},
	{Name => "Daoji Huang",
	    Email => "huan0664@umn.edu",
	    HomePage => "https://daojihuang.me"},
        {Name => "Patricia Klein", 
            Email => "pjklein@tamu.edu", 
            HomePage => "https://patriciajklein.github.io/"},
	{Name => "Adam LaClair",
	    Email => "alaclair@purdue.edu",
	    HomePage=> "https://sites.google.com/view/adamlaclair/home"},
        {Name => "Yuyuan Luo",
            Email => "lyuyuan@mit.edu",
            HomePage=> "https://www.mit.edu/~lyuyuan/"}, 
	{Name => "Joseph McDonough",
	    Email => "mcdo1248@umn.edu",
	    HomePage=> "https://jmcdonough98.github.io/"}
    },
    Headline => "functions for investigating ASM and matrix Schubert varieties",
    PackageExports => {
        "Depth",
        "SimplicialComplexes",
        "SimplicialDecomposability",
        "Posets",
        "MinimalPrimes"
    },
    DebuggingMode => true
)

export{
    
    --MatrixSchubertConstructions.m2
    "isPartialASM",    	       	    --documented ++
    "partialASMToASM",	      	    --documented ++
    "antiDiagInit",    	       	    --documented ++
    "rankTable",    	     	    --documented ++  
    "rotheDiagram",    	       	    --documented ++  
    "augmentedRotheDiagram",	    --documented ++
    "essentialSet",	     	        --documented ++
    "augmentedEssentialSet",        --documented ++
    "schubDetIdeal",	       	    --documented ++
    "fultonGens",    	     	    --documented ++
    "diagLexInitSE",   	      	    --documented ++
    "diagLexInitNW",	    	    --documented ++
    "diagRevLexInit",	     	    --documented ++
    "subwordComplex",	     	    --documented ++
    "entrywiseMinRankTable",	    --documented ++
    "entrywiseMaxRankTable",	    --documented ++
    "schubDecomposition",	        --documented ++
    "permOverASM",                  --documented ++
    "isIntersectionSchubIdeals",    --documented ++
    "isASMIdeal",    	     	    --documented ++	 
    "isASMUnion",    	     	    --documented ++
    "getASM",	     	     	    --documented ++
    "isMinRankTable",	     	    --documented ++
    "rankTableToASM",	     	    --documented ++
    "rankTableFromMatrix",    	    --documented ++
    "schubIntersect",	     	    --documented ++
    "schubAdd",	       	       	    --documented ++
--    "getPermFromASM",	       	    --documented ++
--    "ASM",    	      	      	    -- ??
    
 --permutationMethods.m2   
    "isPerm",	     	     	    --documented ++
    "permToMatrix",    	       	    --documented ++
    "lastDescent",    	      	    --documented ++
    "firstDescent",    	       	    --documented ++
    "permLength",    	     	    --documented ++
    "inverseOf",             	    --documented (check)
    "longestPerm",    	      	    --documented (check)
 --   "getOneReducedWord",    	    -- ??
    "toOneLineNotation",    	    --documented ++
    "composePerms",    	       	    --documented ++
    "isPatternAvoiding",    	    --documented ++
    "isVexillary",    	      	    --documented ++
    "avoidsAllPatterns",	        --documented ++
    "isCartwrightSturmfels",	    --documented ++
    "isCDG",	    	    	    --documented ++
    "rajCode",	      	      	    --documented ++
    "rajIndex",        	       	    --documented ++
    "grothendieckPoly",	       	    -- CHECK DOC
    "schubertPoly",    	       	    -- CHECK DOC
    "doubleSchubertPoly",           -- CHECK DOC
--    "dividedDifference",    	    -- CHECK DOC
--    "PolyType",	       	       	    -- ??
--    "Operator",	       	       	    -- ??
--    "Double",	     	     	    -- ??
    "pipeDreams",    	     	    -- CHECK DOC
    "pipeDreamsNonReduced",    	    -- CHECK DOC
    "netPD",	    	    	    -- CHECK DOC
    "ASMToMonotoneTriangle",        --documented ++
    "MonotoneTriangleToASM",        --documented ++

--MatrixSchubertInvariants.m2    
    "schubReg",    	                --documented ++
    "schubCodim",           	    --documented ++
    "KPolynomialASM",	     	    -- CHECK DOC
    "isSchubCM",    	       	    --documented ++

--ASM_Lists.m2
    "ASMFullList",    	      	    --ADD EX TO DOC
    "ASMRandomList",	    	    --ADD EX TO DOC
    "cohenMacaulayASMsList",	    --ADD EX TO DOC
    "nonCohenMacaulayASMsList",	    --ADD EX TO DOC
    "initialIdealsList"    	        --ADD EX TO DOC
}

--keys used for Schubert/Grothendieck polynomials
protect PolyType
protect Double
protect Operator
protect ASM
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **CODE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
load "./MatrixSchubert/permutationMethods.m2"
load "./MatrixSchubert/MatrixSchubertConstructions.m2"
load "./MatrixSchubert/MatrixSchubertInvariants.m2"
load "./MatrixSchubert/ASM_Lists.m2"
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **DOCUMENTATION** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
beginDocumentation ()    
load "./MatrixSchubert/permutationMethodsDOC.m2"
load "./MatrixSchubert/MatrixSchubertConstructionsDOC.m2"
load "./MatrixSchubert/MatrixSchubertInvariantsDOC.m2"
load "./MatrixSchubert/ASM_ListsDOC.m2"
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **TESTS** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
load "./MatrixSchubert/MatrixSchubertTests.m2"

end---------------------------------------------------------------------------     

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **SCRATCH SPACE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------


---------------------
--Ayah's sandbox
---------------------


------------------------------------
--Development Section
------------------------------------

restart
uninstallPackage "MatrixSchubert"
restart
installPackage "MatrixSchubert"
restart
needsPackage "MatrixSchubert"
elapsedTime check "MatrixSchubert"
viewHelp "MatrixSchubert"
