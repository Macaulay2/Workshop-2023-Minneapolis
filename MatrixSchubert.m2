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

--TODO: organize the exports
export{
    "isPartialASM",
    "partialASMToASM",
    "fultonGens",
    "schubertDetIdeal",
    "diagLexInit",
    "antiDiagInit",
    "rankMatrix",
    "essentialBoxes",
    "subwordComplex",
    "grothendieckPoly",
    "rotheDiagram",
    "permToMatrix",
    "toOneLineNotation",
    "composePerms",
    "isPerm",
    "longestPerm",
    "dividedDifference",
    "schubertPoly",
    "doubleSchubertPoly",
    "permLength",
    "augmentedRotheDiagram",
    "isPatternAvoiding",
    "isVexillary",
    "isCDG",
    "inverseOf",
    "isCartwrightSturmfels",
    "schubertDecomposition",
    "isIntersectionSchubIdeals",
    "rajCode",
    "rajIndex",
    "isMinRankTable",
    "Double",
    "Operator",
    "PolyType",
    "rankTableToASM",
    "rankTableFromMatrix",
    "schubertCodim",
    "matrixSchubertRegADI",
    "matrixSchubertReg",
    "isASMIdeal",
    "ASM",
    "getASM",
    "isSchubertCM",
    "ASMFullList",
    "cohenMacaulayASMsList",
    "nonCohenMacaulayASMsList",
    "idealsList",
    "lastDescent",
    "firstDescent"
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

    apply( ran(toList(1..n)),w->assert(matrixSchubertRegADI(w) == matrixSchubertReg(w)));
        
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

------------------------------------
--Development Section
------------------------------------

restart
uninstallPackage "MatrixSchubert"
restart
installPackage "MatrixSchubert"
restart
debug needsPackage "MatrixSchubert"
elapsedTime check "MatrixSchubert"
viewHelp "MatrixSchubert"

