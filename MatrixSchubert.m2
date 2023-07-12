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
	    Email => "",
	    HomePage => "daojihuang.me"},
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
	    HomePage=> " "}
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
    "rankMatrix",    	     	    --documented ++  
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
    "entrywiseMinRankTable",	    --documented (check)
    "entrywiseMaxRankTable",	    --documented (check)
    "schubDecomposition",	    --documented ++
    "permOverASM",                  --documented ++
    "isIntersectionSchubIdeals",    -- ADD EX TO DOC
    "isASMIdeal",    	     	    -- ADD EX TO DOC	 
    "isASMUnion",
    "getASM",	     	     	    -- ADD EX TO DOC
    "isMinRankTable",	     	    --documented ++
    "rankTableToASM",	     	    --documented ++
    "rankTableFromMatrix",    	    --documented ++
    "schubIntersect",	     	    -- documented (check)
    "schubAdd",	       	       	    -- documented (check)
    "getPermFromASM",	       	    --documented (check)
    "ASM",    	      	      	    -- ??
    
 --permutationMethods.m2   
    "isPerm",	     	     	    --documented ++
    "permToMatrix",    	       	    --documented ++
    "lastDescent",    	      	    --documented ++
    "firstDescent",    	       	    --documented ++
    "permLength",    	     	    --documented ++
    "inverseOf",             	    -- ??
    "longestPerm",    	      	    -- ??
    "getOneReducedWord",    	    -- ??
    "toOneLineNotation",    	    --documented ++
    "composePerms",    	       	    --documented ++
    "isPatternAvoiding",    	    --documented ++
    "isVexillary",    	      	    --documented ++
    "avoidsAllPatterns",	       	    -- CHECK DOC
    "isCartwrightSturmfels",	    -- CHECK DOC
    "isCDG",	    	    	    -- CHECK DOC
    "rajCode",	      	      	    --documented ++
    "rajIndex",        	       	    --documented ++
    "grothendieckPoly",	       	    -- CHECK DOC
    "schubertPoly",    	       	    -- CHECK DOC
    "doubleSchubertPoly",           -- CHECK DOC
    "dividedDifference",    	    -- CHECK DOC
    "PolyType",	       	       	    -- ??
    "Operator",	       	       	    -- ??
    "Double",	     	     	    -- ??

--MatrixSchubertInvariants.m2    
    "schubReg",    	            --documented ++
    "schubCodim",           	    --documented ++
    "KPolynomialASM",	     	    -- ??

--Lists.m2
    "ASMFullList",    	      	    --ADD EX TO DOC
    "ASMRandomList",	    	    --ADD EX TO DOC
    "cohenMacaulayASMsList",	    --ADD EX TO DOC
    "nonCohenMacaulayASMsList",	    --ADD EX TO DOC
    "initialIdealsList"    	    --ADD EX TO DOC
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
debug needsPackage "MatrixSchubert"

w = {2,1,4,3}
I = schubertDetIdeal w;
R = ring I;
kk = coefficientRing R;
apply(#w, i-> toList insert(i,1,(#w-1):0))
degs = splice apply(oo, i-> #w:i)
Q = kk[R_*, Degrees => degs]
G = numerator hilbertSeries sub(I, Q)
factor G
S =  ring G
F = map(S,S, {T_0=>1-T_0, T_1 => 1-T_1, T_2 => 1-T_2})
F G
factor oo
schubertPoly {2,1,4,3}

M = matrix{{0,0,1},{1,0,-1}}
KPolynomialASM M
-----------------------------------------
--Adam's Testing for matrixSchubertReg --
-----------------------------------------


Tester (n) -> (

    S = apply(permutations(n),S->apply(S,i->i+1));
    assert(apply(S,i->matrixSchubertRegADI(i))==apply(S,i->matrixSchubertRegWPS(i)))
);

TesterPerm = (w) -> (
    assert(matrixSchubertRegADI(w)==matrixSchubertReg(w));
);

NTests = 50;
NMax = 5;
for i from 1 to NTests do (
    w = random toList (1..(random(1,NMax)));
    TesterPerm(w);
    << "Finished test #" << i << << ", perm size: " << length(w) << endl;
);

randomPerm = (Len) -> (
    return random toList (1..random(1,Len));
);

TesterPermTime = (w) -> (
    elapsedTime matrixSchubertRegWPS w;
--    elapsedTime matrixSchubertRegADI w;
    << endl;
);

Tester = (n) -> (

    --apply( random(toList(1..n)),w->assert(schubReg(w) == schubReg(permToMatrix(w))));
    w =  random(toList(1..n));
    assert(schubReg(w) == schubReg(permToMatrix(w)));
);

randomTester = (lenPerm) -> (

    w = random (toList(1..lenPerm));
    
    << w << endl;

    elapsedTime (rADI = matrixSchubertRegADI(w));
    elapsedTime (rRaj = matrixSchubertReg(w));

    if (rADI != rRaj) then (
	<< "BAD: " << w << endl;
	exit(1);
    );
);

apply(1..10,i->Tester(i));



------------------------------------
--Adam RajCode Testing--
------------------------------------

rajTest = (w) -> (
    if(rajCode(w) != rajCodeRec(w)) then (
    	print(w);
	exit ;	
    );
);

randomTest = (n) -> (
  
    w= random toList(1..n);
    rajTest(w)
--    print elapsedTime rajCode(w);
    print elapsedTime rajCodeRec(w);
    << endl;
);

for i from 1 to 4 do (
    for w in permutations(toList (1..i)) do rajTest(w);
);

setRandomSeed 50;

for i from 1 to 10 do (
     i = random toList (1..7);
     << i << endl;
     print elapsedTime matrixSchubertRegADI(i);
     << endl << endl;
);

for l in permutations(toList(1..4)) do (
    --print permToMatrix l;
      if (schubReg permToMatrix l != 0) then print l;
      );

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



