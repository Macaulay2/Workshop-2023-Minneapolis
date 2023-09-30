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
	    Email => "alaclair@purdue.edu ",
	    HomePage=> "https://sites.google.com/view/adamlaclair/home"},
        {Name => "Yuyuan Luo",
            Email => "lyuyuan@mit.edu",
            HomePage=> "https://www.mit.edu/~lyuyuan/"}, 
	{Name => "Joseph McDonough",
	    Email => "mcdo1248@umn.edu",
	    HomePage=> ""}
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
    "essentialSet",	     	    --documented ++
    "augmentedEssentialSet",        --documented ++
    "schubDetIdeal",	       	    --documented ++
    "fultonGens",    	     	    --documented ++
    "diagLexInitSE",   	      	    --documented ++
    "diagLexInitNW",	    	    --documented ++
    "diagRevLexInit",	     	    --documented ++
    "subwordComplex",	     	    --documented ++
    "entrywiseMinRankTable",	    --documented ++
    "entrywiseMaxRankTable",	    --documented ++
    "schubDecomposition",	    --documented ++
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
    "getPermFromASM",	       	    --documented ++
    "ASM",    	      	      	    -- ??
    
 --permutationMethods.m2   
    "isPerm",	     	     	    --documented ++
    "permToMatrix",    	       	    --documented ++
    "lastDescent",    	      	    --documented ++
    "firstDescent",    	       	    --documented ++
    "permLength",    	     	    --documented ++
    "inverseOf",             	    --documented (check)
    "longestPerm",    	      	    --documented (check)
    "getOneReducedWord",    	    -- ??
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
    "dividedDifference",    	    -- CHECK DOC
    "PolyType",	       	       	    -- ??
    "Operator",	       	       	    -- ??
    "Double",	     	     	    -- ??
    "pipeDreams",    	     	    -- ??
    "pipeDreamsNonReduced",    	    -- ??
    "netPD",	    	    	    -- ??
    "ASMToMonotoneTriangle",        --documented ++
    "MonotoneTriangleToASM",        --documented ++

--MatrixSchubertInvariants.m2    
    "schubReg",    	                --documented ++
    "schubCodim",           	    --documented ++
    "KPolynomialASM",	     	    -- ??
    "isSchubertCM",    	       	    --documented ++

--Lists.m2
    "ASMFullList",    	      	    --ADD EX TO DOC
    "ASMRandomList",	    	    --ADD EX TO DOC
    "cohenMacaulayASMsList",	    --ADD EX TO DOC
    "nonCohenMacaulayASMsList",	    --ADD EX TO DOC
    "initialIdealsList"    	        --ADD EX TO DOC
}
    
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **CODE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
load "./MatrixSchubert/Code/helpers.m2"
load "./MatrixSchubert/Code/permutationMethods.m2"
load "./MatrixSchubert/Code/MatrixSchubertConstructions.m2"
load "./MatrixSchubert/Code/MatrixSchubertInvariants.m2"
load "./MatrixSchubert/Code/lists.m2"
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **DOCUMENTATION** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
beginDocumentation ()    
load "./MatrixSchubert/Documentation/permutationMethodsDOC.m2"
load "./MatrixSchubert/Documentation/MatrixSchubertConstructionsDOC.m2"
load "./MatrixSchubert/Documentation/MatrixSchubertInvariantsDOC.m2"
load "./MatrixSchubert/Documentation/listsDOC.m2"
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **TESTS** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
load "./MatrixSchubert/Tests/MatrixSchubertTests.m2"
load "./MatrixSchubert/Tests/MatrixSchubertTestsIdentity.m2"
load "./MatrixSchubert/Tests/MatrixSchubertTestsInterestingExamples.m2"

end---------------------------------------------------------------------------     

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **SCRATCH SPACE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------


---------------------
--Ayah's sandbox
---------------------

restart
uninstallPackage "MatrixSchubert"
restart
installPackage "MatrixSchubert"
restart
needsPackage "MatrixSchubert"
initialIdealsList 3

    n = 3
    if (n < 3 or n > 6) then error("There is no available list for this n.");
    filename = concatenate("./MatrixSchubert/ASMData/antiDiagIniIdeal/ideals", toString n, ".txt");
    z = getSymbol "z";
    S = QQ(monoid[z_(1,1)..z_(n,n)]);
    listOfIdeals = apply(lines get filename, i -> (use S; value i))
    listOfIdeals

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

