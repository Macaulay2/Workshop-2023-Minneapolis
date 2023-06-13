newPackage(
    "FInvariants",
    Version => "1.1",
    Date => "June 8, 2023",
    Authors => {{Name => "Anna Brosowsky", Email => "annabro@umich.edu"},
	{Name => "Havi Ellers", Email => "ellers@umich.edu"},
	{Name => "Moty Katzman", Email => "M.Katzman@shef.ac.uk"},
	{Name => "Jiamin Li", Email => "jli283@uic.edu"},
	{Name => "Swaraj Pande", Email => "swarajsp@umich.edu"},
	{Name => "Abraham Pascoe", Email => "abrahampascoe@ku.edu"},
	{Name => "Austyn Simpson", Email => "austyn@umich.edu"},
	{Name => "Pedro Teixeira", Email => "pteixeir@knox.edu"}
        },
    Headline => "A package for calculations of Frobenius invariants",
    DebuggingMode => true,
    Reload => true,
    AuxiliaryFiles => true,
    PackageExports => { "TestIdeals" },
    PackageImports => { "Saturation", "PrimaryDecomposition" }
)

importFrom( "ReesAlgebra", { "Tries" } )

export {
    -- in this file
    "HSL", 
    "FrobeniusonExt", 
    "HSLNumber",
    -- in F-modules
    "FModule",
    "GeneratingMorphism",
    "FF",
    "generatingMorphism",
    "makeFModule",
    "localCohomology",
    "root",
    "rootMorphism",
    "cohomDim",
    "randomGeneratingMorphism",
    "isFilterRegElement",
    "isFilterRegSeq",
    "randomFilterRegSeq",
    "limitClosure",
    "lowerLimit",
    "lyubeznikNumber",
    "lyubeznikTable",
    "MaxDegree"
}

load "./FInvariants/F-Modules.m2"

HSL = method()

HSL ( Module, Matrix ) := ( A, u ) -> 
(
    A1 := matrix entries presentation A;
    u1 := matrix entries u;
    e := 0;
    M := gens target u1;
    local newM;
    while true do 
    (
	newM = matrix entries frobeniusRoot( 1, u1*M );
	if image( newM ) + image( A1 ) == image( M ) + image( A1 ) then return e;
	e = e + 1;
	M = newM;
    )
)

HSL ( Matrix, Matrix ) := ( A, u ) -> HSL( coker A, u )

FrobeniusonExt = method()

-- This function takes an ideal in a polynomial ring and a local cohomology index and returns the frobenius map on the Matlis dual of local Cohomology of the quotient module.
FrobeniusonExt ( Ideal, ZZ ) := ( I, i ) -> 
(
    R := ring I;
    d := dim R;   
    M2 := R^1/I;
    j := d - i;
    M1 := R^1/(frobenius I);
    phi := map( M2, M1, 1 );
    Ext^j( phi, R^1 )
)

HSLNumber = ( I, i ) -> 
(
    R := ring I;
    d := dim R;
    j := d - i;
    u := FrobeniusonExt( I, i );
    A := Ext^j( R^1/I, R^1 );
    HSL( A, u )
)

beginDocumentation()

----------------------------------------------------------------------
-- Section for Tests ------------------------------------------------
---------------------------------------------------------------------

-- To create a test, use the format TEST <string>, where string is a string
-- containing valid M2 code for a test. Note that /// is alternate string syntax
-- to ", and is probably a better idea since then safe to include quotes

-- To run these tests, simply open m2 (probably from this folder), and do:
-- needsPackage "FInvariants"
-- check FInvariants

load "./FInvariants/F-Modules_test.m2"

--Test example, should output 2
TEST ///
--needsPackage "FInvariants"
R = ZZ/2[x1,x2,x3,x4,x5];
I = ideal(x2^2+x1*x3,x1*x2*x4^2+x3^3*x5,x1^2*x4^2+x2*x3^2*x5);
M2 = R^1/I;
A = Ext^2( M2, R^1 );
M1 = R^1/( frobenius I );
phi = map( M2, M1, 1 );
u = Ext^2( phi, R^1 );
B = presentation A;
assert( HSL( B, u ) == 2 ) -- this is the thing we want to be true
///

--Test example, should output 1
TEST ///
-- needsPackage "FInvariants"
R = ZZ/11[x,y,z];
I = ideal( x^5+y^5+z^5 );
M2 = R^1/I;
A = Ext^1( M2, R^1 );
M1 = R^1/( frobenius I );
phi = map( M2, M1, 1 );
u = Ext^1( phi, R^1 );
B = presentation A;
assert( HSL( B, u ) == 1 )
///

-- to see what failure looks like, uncomment this fake function
-- TEST /// assert(1==2); ///

end--
