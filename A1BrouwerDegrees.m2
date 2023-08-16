--A1BrouwerDegrees.m2
newPackage(
    "A1BrouwerDegrees",
    Version=>"1.0",
    Date=>"June 5, 2023",
    Authors=>{
        {Name=>"Nikita Borisov",
	 Email=>"nborisov@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/nikita-borisov"},
        {Name=>"Thomas Brazelton",
	 Email=>"tbraz@math.upenn.edu",
	 HomePage=>"https://www2.math.upenn.edu/~tbraz/"},
        {Name=>"Frenly Espino",
	 Email=>"frenly@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/frenly-espino"},
         {Name=>"Tom Hagedorn",
	 Email=>"hagedorn@tcnj.edu",
	 HomePage=>"https://hagedorn.pages.tcnj.edu/"},
        {Name=>"Zhaobo Han",
	 Email=>"zbtomhan@sas.upenn.edu",
	 HomePage=>"https://www.linkedin.com/in/zhaobo-han-77b1301a2/"},
     	{Name=>"Jordy Lopez Garcia",
	 Email=>"jordy.lopez@tamu.edu",
	 HomePage=>"https://jordylopez27.github.io/"},
        {Name=>"Joel Louwsma",
	 Email=>"jlouwsma@niagara.edu",
	 HomePage=>"https://www.joellouwsma.com/"},
        {Name=>"Andrew Tawfeek",
	 Email=>"atawfeek@uw.edu",
	 HomePage=>"https://www.atawfeek.com/"},
        {Name=>"Wern Juin Gabriel Ong",
	 Email=>"gong@bowdoin.edu",
	 HomePage=>"https://wgabrielong.github.io/"}
	},
    Headline=>"TODO",
    PackageImports=>{
	"Parametrization",
	"RealRoots",
	"RationalPoints2"
	},
    PackageExports=>{},
    DebuggingMode=>true
    )


export{
    -- ArithmeticMethods.m2
    "squarefreePart",
    "legendreBoolean",
    "localAlgebraBasis",
    "primeFactors",
    "PadicValuation",
    "squareSymbol",

    --MatrixMethods.m2
    "congruenceDiagonalize",
    
    --GrothendieckWittClasses.m2    
    "GrothendieckWittClass",
    "baseField",
    "gwClass",
    "gwAdd",
    "gwMultiply",
    
    --BuildingForms.m2
    "diagonalClass",
    "hyperbolicForm",
    "PfisterForm",
    
    --SimplifiedRepresentatives.m2
    "diagonalForm",
    "diagonalEntries",
    "integralDiagonalRep",
    
    --HilbertSymbols.m2
    "HilbertSymbol",
    
    --GWInvariants.m2
    "signature",
    "integralDiscriminant",
    "relevantPrimes",
    "HasseWittInvariant",

    --LocalGlobalDegrees.m2
    "globalA1Degree",
    "localA1Degree",
    
    --IsomorphismOfForms.m2
    "isIsomorphicFormQ",
    "gwIsomorphic",
    
    --Isotropy.m2
    --"isIsotropicQp",
    "isIsotropic",
    "isAnisotropic",

    --AnisotropicDimension.m2
    "isHyperbolicQp",
    "anisotropicDimensionQp",
    "anisotropicDimensionQQ",
    "anisotropicDimension",
    "WittIndex",
    
    --Decomposition.m2
    "sumDecomposition",
    "sumDecompositionString"
    }

-- Basic arithmetic, p-adic, and commutative algebra operations we will use
load "./A1-Brouwer/Code/ArithmeticMethods.m2"

-- Basic manipulations of matrices we will use
load "./A1-Brouwer/Code/MatrixMethods.m2"

-- Establishing the GrothendieckWittClass type and some basic manipulations
load "./A1-Brouwer/Code/GrothendieckWittClasses.m2"

-- For building new symmetric bilinear forms
load "./A1-Brouwer/Code/BuildingForms.m2"

-- For providing simplified representatives of symmetric bilinear forms
load "./A1-Brouwer/Code/SimplifiedRepresentatives.m2"

-- For Hilbert symbols over p-adic numbers
load "./A1-Brouwer/Code/HilbertSymbols.m2"

-- Invariants of symmetric bilinear forms
load "./A1-Brouwer/Code/GWInvariants.m2"
    
-- Local and global A1-brouwer degrees
load "./A1-Brouwer/Code/LocalGlobalDegrees.m2"

-- Checking if forms are isomorphic
load "./A1-Brouwer/Code/IsomorphismOfForms.m2"

-- For verifying (an)isotropy
load "./A1-Brouwer/Code/Isotropy.m2"

-- Anisotropic dimension
load "./A1-Brouwer/Code/AnisotropicDimension.m2"

-- Finally, decomposing forms
load "./A1-Brouwer/Code/Decomposition.m2"

load "./A1-Brouwer/MoreMatrixMethods.m2"
load "./A1-Brouwer/NondegeneratePartDiagonal.m2"






----------------------------
----------------------------
-- DOCUMENTATION
----------------------------
----------------------------

beginDocumentation()

document{
    Key => A1BrouwerDegrees,
    Headline => "A package for running A1-Brouwer degree computations in Macaulay2",
    }

undocumented {
    }

load "./A1-Brouwer/Documentation/ArithmeticMethodsDoc.m2"

load "./A1-Brouwer/Documentation/MatrixMethodsDoc.m2"

load "./A1-Brouwer/Documentation/GrothendieckWittClassesDoc.m2"

load "./A1-Brouwer/Documentation/BuildingFormsDoc.m2"

load "./A1-Brouwer/Documentation/SimplifiedRepresentativesDoc.m2"

load "./A1-Brouwer/Documentation/HilbertSymbolsDoc.m2"

load "./A1-Brouwer/Documentation/GWInvariantsDoc.m2"

load "./A1-Brouwer/Documentation/LocalGlobalDegreesDoc.m2"

load "./A1-Brouwer/Documentation/IsomorphismOfFormsDoc.m2"

load "./A1-Brouwer/Documentation/IsotropyDoc.m2"

load "./A1-Brouwer/Documentation/AnisotropicDimensionDoc.m2"

load "./A1-Brouwer/Documentation/DecompositionDoc.m2"












----------------------------
----------------------------
-- Testing
----------------------------
----------------------------

-- For debugging: remember Macaulay2 starts counting the first test as Test 0

       
------------------
-- TESTING
------------------

-- Diagonal form testing
TEST ///
print("diagonal form testing");
M1=matrix(RR, {{0, 1}, {1, 0}});
G1=gwClass(M1);
M2=diagonalForm(G1);
assert(M2.matrix===matrix(RR, {{1, 0}, {0, -1}}));
///

TEST ///
M3=matrix(CC, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
assert(M4.matrix===matrix(CC, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
-- Ensure the cache is populated
-- assert(G2.cache.?diagonalForm)
///


TEST ///
M3=matrix(QQ, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
assert(M4.matrix===matrix(QQ,{{1, 0, 0}, {0, -2, 0}, {0, 0, 1/2}}));
///

-- gwTypeTest.m2
TEST ///
M = matrix(QQ,{{1,0},{0,1}});
N = matrix(QQ, {{1, 2}, {2, 5}})
beta = gwClass(M);
gamma = gwClass(N);
assert(baseField(beta) === QQ)
assert(beta.matrix === M)
--Operations within GW-classes
A = gwAdd(beta, gamma);
B = gwMultiply(beta, gamma);
assert(A.matrix === matrix(QQ, {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 2}, {0, 0, 2, 5}}));
assert(B.matrix === matrix(QQ, {{1, 2, 0, 0}, {2, 5, 0, 0}, {0, 0, 1, 2}, {0, 0, 2, 5}}));
---non well-defined GW-classes
M'=matrix(ZZ, {{1, 0}, {0, 1}});
N'=matrix(QQ, {{1, 1}, {1, 1}});
theta=gwClass(M');
sigma=gwClass(N');
assert(isWellDefined(theta) === false);
assert(isWellDefined(sigma) === false);
///


-- degreeTesting.m2
TEST ///
T1 = QQ[x]
f = {x^2}
beta = globalA1Degree(f)
gamma = gwClass(matrix(QQ,{{0,1},{1,0}}))
assert(gwIsomorphic(beta,gamma))
///


TEST ///
T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
f1GD = globalA1Degree(f1);
f1GDmat = f1GD.matrix;
assert((WittDecomp(f1GDmat))==(3,matrix(QQ,{{}})));
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
f1LDq= localA1Degree(f1,q);
f1LDr= localA1Degree(f1,r);
f1LDsum = gwAdd(f1LDq, f1LDr);
assert(gwIsomorphic(f1LDsum, f1GD));



T2 = QQ[w];
f2 = {w^4 + w^3 - w^2 - w};
f2GD= globalA1Degree(f2);
f2GDmat = f2GD.matrix;
assert(WittDecomp(f2GDmat)==(2,matrix(QQ,{{}})));

p=ideal {w+1};
f2LDp = localA1Degree(f2, p);
f2LDpmat = f2LDp.matrix;
assert(WittDecomp(f2LDpmat)==(1,matrix(QQ,{{}})));
s=ideal{w-1};
f2LDs = localA1Degree(f2, s);
t=ideal{w};
f2LDt = localA1Degree(f2, t);
f2LDsum = gwAdd(gwAdd(f2LDp, f2LDs),f2LDt);
assert(gwIsomorphic(f2LDsum, f2GD));
///
    


end