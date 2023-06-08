newPackage(
         "RandomIdealsPlus",
         Version => "0.1",
         Date => "",
	 AuxiliaryFiles => false,
         Authors => {{ Name => "", Email => "", HomePage => ""}},
	 Headline => "Creating random ideals for experiments",
	 PackageExports => {"DGAlgebras", "RandomIdeals"},
         DebuggingMode => false )
export {
    "mixedIdeal",
    "mixedIdealList"
}

-* Code section *-
  

        
mixedIdeal = method()
mixedIdeal (Ring, List, List) := Ideal => (S, BList, MList)-> (
 randomPureBinomialIdeal(BList, S)+ randomMonomialIdeal(MList,S)
 )   	


mixedIdealList = method()
mixedIdealList(Ring, ZZ, List, List) := List => (A, Num, Bdegs, Mdegs) -> (
    for i from 1 to Num list (
	trim(randomBinomialIdeal(Bdegs,A)+randomMonomialIdeal(Mdegs,A))
	)
    )


-* Documentation section for mixedIdeal *-
beginDocumentation()

doc ///
Key
 RandomIdealsPlus
Headline
 Functions to produce Random ideals. 
Description 
 Text
  Collection of functions  to create random ideals (or) a list of random ideals
  with some constraints imposed.
///


doc ///
Key
  mixedIdeal
Headline
 Computes an ideal which is the sum of a binomial ideal and a monomial ideal with the generators of specified degrees. 
Usage
  I = mixedIdeal(R,L1,L2)   
Inputs
 R: Ring
 L1: List
 L2: List
Outputs
 I: Ideal
Description
  Text
    Mixed ideal(in this context) is the sum of a monomial ideal and a binomial ideal.
  Example
    R = ZZ/101[a,b,c]
    Blist = {3,4}
    Mlist = {2,3}
    I = mixedIdeal(R, Blist, Mlist)
///


doc ///
Key
  mixedIdealList
Headline
 Computes a list of ideals which is the sum of a binomial ideal and a monomial ideal with the generators of specified degrees.
Usage
  L  = mixedIdealList(R,i,L1,L2)   
Inputs
 R: Ring
 i: ZZ
 L1: List
 L2: List
Outputs
 L: List
Description
  Text
    The function computes a list of ideals which is the sum of a binomial ideal and a monomial ideal with the 
    generators of specified degrees.
  Example
    R = ZZ/101[a,b,c]
    Blist = {3,4}
    Mlist = {2,3}
    I = mixedIdealList(R,3,Blist, Mlist)
///





-* Test section *-
TEST /// -* [insert short title for this test] *-
-- test code and assertions here
-- may have as many TEST sections as needed
///

end--

-* Development section *-
restart
loadPackage("RandomIdealsPlus", Reload => true)
debug needsPackage "RandomIdealsPlus"
check "RandomIdealsPlus"

uninstallPackage "RandomIdealsPlus"
restart
installPackage "RandomIdealsPlus"
viewHelp "RandomIdealsPlus"

	"randomMonomial",
	"randomArtinIdeal",
	"randomArtinDegreeIdeal",
	"randomArtinMonomialIdeal",
	"randomArtinDegreeMonomialIdeal",
	"randomArtinRing",
	"randomArtinDegreeRing",
	"randomArtinMonomialRing",
	"randomArtinDegreeMonomialRing",
	"random toric ideal"


