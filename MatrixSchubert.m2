newPackage(
        "MatrixSchubert",
	AuxiliaryFiles => true,
        Version => "0.1", 
        Date => "",
        Authors => {
	    {Name => "Ayah Almousa", 
                Email => "almou007@umn.edu", 
                HomePage => "http://sites.google.com/view/ayah-almousa"},
            {Name => "Patricia Klein", 
                Email => "pjklein@tamu.edu", 
                HomePage => " "},
	    {Name => "Yuyuan Luo",
		Email => "",
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
    "isPartialASM",
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
    "composePerms",
    "isPerm",
    "schubertPoly",
    "doubleSchubertPoly",
    "entrywiseMinRankTable",
    "entrywiseMaxRankTable",
    "permLength",
    "augmentedRotheDiagram",
    "isPatternAvoiding",
    "isVexillary",
    "schubertDecomposition",
    "isIntersectionSchubIdeals",
    "rajCode",
    "rajIndex",
    "isMinRankTable",
    "Double",
    "rankTableToASM",
    "schubertReg",
    "rankTableFromMatrix",
    "schubertCodim",
    "matrixSchubertRegADI",
    "matrixSchubertRegWPS",
    "lengthIncrSubset",
    "isASMIdeal",
    "ASM",
    "getASM",
    "isSchubertCM",
    "ASMFullList"
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

M = matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}}
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


A = matrix{{0,-1,0,1,1},{1,-1,1,-1,1},{0,1,1,0,-1},{1,1,-1,1,-1},{-1,1,0,0,1}}
n = numrows A
m = numcols A
rowCheck = new MutableList
colCheck = new MutableList
for i from 0 to n-1 do(
    partialSums = for i from 0 to m-1 list(sum(delete(0, flatten entries A_{i})));
    rowCheck#(#rowCheck) = (((unique sort partialSums) == {0,1}) or ((unique sort partialSums) == {0}) or ((unique sort partialSums) == {1}));
    );
for i from 0 to m-1 do(
    partialSums = for i from 0 to n-1 list(sum(delete(0, flatten entries A^{i})));
    rowCheck#(#colCheck) = (((unique sort partialSums) == {0,1}) or ((unique sort partialSums) == {0}) or ((unique sort partialSums) == {1}));
    );
(toList rowCheck == toList(#rowCheck:true)) and (toList colCheck == toList(#colCheck:true))
)



-----------------------------------------
--Adam's Testing for matrixSchubertReg --
-----------------------------------------


Tester (n) -> (

    S = apply(permutations(n),S->apply(S,i->i+1));
    assert(apply(S,i->matrixSchubertRegADI(i))==apply(S,i->matrixSchubertRegWPS(i)))
);

TesterPerm = (w) -> (
    assert(matrixSchubertRegADI(w)==matrixSchubertRegWPS(w));
);

NTests = 50;
NMax = 8;
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

apply(1..10,i->Tester(i));

restart
needsPackage "MatrixSchubert"
I = schubertDetIdeal {2,1,6,3,5,4}
isASMIdeal I
I.cache.ASM
getASM(I)

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


