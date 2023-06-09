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
	    HomePage => " "},
        {Name => "Patricia Klein", 
            Email => "pjklein@tamu.edu", 
            HomePage => " "},
	{Name => "Adam LaClair",
	    Email => "",
	    HomePage=> " "},
        {Name => "Yuyuan Luo",
            Email => "",
            HomePage=> " "},
	{Name => "Joe McDonough",
	    Email => " ",
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
    "schubertPoly",
    "doubleSchubertPoly",
    "permLength",
    "augmentedRotheDiagram",
    "isPatternAvoiding",
    "isVexillary",
    "isCDG",
    "isCartwrightSturmfels",
    "schubertDecomposition",
    "isIntersectionSchubIdeals",
    "rajCode",
    "rajIndex",
    "isMinRankTable",
    "Double",
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
    "lastDescent",
    "firstDescent",
    "indexOfVariable" -- todo, figure out a way to test this without exporting it
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

M = matrix{{0,0,1,0,0},{0,1,-1,1,0}}
fultonGens M
schubertDetIdeal M

w = {3,2,5,1,4}
schubertDetIdeal w

schubertDetIdeal {1,2,3}


permToMatrix w
isPartialASM M
rotheDiagram M
essentialBoxes M
grothendieckPoly M
netList fultonGens M
subwordComplex M
betti res diagLexInit M
betti res antiDiagInit M


w = {4,1,3,5,2}
A = (permToMatrix w)_{0,1,2}

if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
zMatrix = genMat(numrows A, numcols A); --generic matrix
rankMat = rankMatrix A; --rank matrix for A
essBoxes = essentialBoxes A;
if essBoxes == {} then (
    R = ring zMatrix;
    return ideal(0_R)
    );
zBoxes = apply(essBoxes, i-> flatten table(i_0,i_1, (j,k)->(j+1,k+1))); --smaller matrix indices for each essential box
ranks = apply(essBoxes, i-> rankMat_(i_0-1,i_1-1)); --ranks for each essential box
antiDiagGens = new MutableList;
for box in essBoxes do (
    pos = position(essBoxes, i-> i==box);
    boxSubmatrix = zMatrix^{0..(box_0-1)}_{0..(box_1-1)};
    for x in subsets(numrows boxSubmatrix,ranks_pos+1) do (
    	for y in subsets(numcols boxSubmatrix,ranks_pos+1) do(
	    indicesList = apply(pack(2,mingle(x,reverse y)),i->toSequence i);
	    antiDiagGens#(#antiDiagGens) =  product(apply(indicesList,i->boxSubmatrix_i));
	    );
	);
    );
ideal(unique flatten toList antiDiagGens)





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
needsPackage "MatrixSchubert"
elapsedTime check "MatrixSchubert"
viewHelp "MatrixSchubert"
