newPackage(
    "MatrixFactorizations",
    AuxiliaryFiles => true,
    Version => "0.1", 
    Date => "",
    Authors => {
        {Name => "David Favero", 
            Email => "favero@umn.edu", 
            HomePage => ""},
	{Name=> "Sasha Pevzner",
	    Email => "pevzn002@umn.edu",
	    HomePage => ""},
	{Name => "Timothy Tribone",
	    Email => "tim.tribone@utah.edu",
	    HomePage => ""},
        {Name => "Keller VandeBogert", 
            Email => "kvandebo@nd.edu", 
            HomePage => ""}
    },
    Headline => "a package for computing with length d matrix factorizations",
    PackageExports => {
        "Complexes",
	"CompleteIntersectionResolutions",
	"TensorComplexes"
    },
    DebuggingMode => true
)

export{
    
    --ZZdFactorizations.m2
    "ZZdFactorization",
    "ZZdFactorizationMap",
    "period",
    "ZZdfactorization",
    "lineOnTop",
    "Fold",
    "component",
    "isCommutativeCached",
    "isZZdFactorizationMorphism",
    "dHom",
    "trans",
    "randomFactorizationMap",
    "tensorCommutativity",
    "isNullHomotopyOf",
    "isNullHomotopic",
    "nullHomotopy",
    "isShortExactSequence",
    
    --ChernCharacter.m2
    "chernCharacter",
    "boundaryBulk",
    
    --CIResCompatibility.m2
    "higherHomotopyFactorization",
    "toBranchedCover",
    "mooreMF",
    "rk1MCM2gen",
    "rk1MCM3gen",
    "classicalAdjoint",
    
    --Compositions.m2
    "posComps",
    "compHash",
    
    --eisenbud-schreyer-examples.m2
    "quadraticMF",
    "ulrichFromMF",
    "polynomial", --This one should be in main file, but was having an issue
    
    --functionsMF_new.m2
    "shift",
    "directSumMF",
    "isdFactorization",
    "tailMF",
    "tensorMF",
    "dTensor",
    "koszulMF",
    "koszulMFf",
    "eulerMF",
    "freeMods",
    "adjoinRoot",
    "trivialFactorization",
    "toComplex",
    "HHMF",
    "ecMF",
    
    --KoszulFactorization.m2
    --don't think this file is needed?
    --tensorMF.m2 seems to also not be needed
}
    
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **CODE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
load "./Matrix-Factorizations/ZZdFactorizations.m2"
load "./Matrix-Factorizations/chernCharacter.m2"
load "./Matrix-Factorizations/CIResCompatibility.m2"
load "./Matrix-Factorizations/Compositions.m2"
load "./Matrix-Factorizations/eisenbud-schreyer-examples.m2"
load "./Matrix-Factorizations/functionsMF_new.m2"
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **DOCUMENTATION** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
beginDocumentation ()    
load "./Matrix-Factorizations/ZZdFactorizationsDOC.m2"
------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **TESTS** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------

end---------------------------------------------------------------------------     

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **SCRATCH SPACE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------


------------------------------------
--Development Section
------------------------------------

restart
uninstallPackage "MatrixFactorizations"
restart
installPackage "MatrixFactorizations"
restart
needsPackage "MatrixFactorizations"
elapsedTime check "MatrixFactorizations"
viewHelp "MatrixFactorizations"
