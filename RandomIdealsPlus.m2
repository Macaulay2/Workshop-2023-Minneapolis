newPackage(
         "RandomIdealsPlus",
         Version => "0.1",
         Date => "",
	 AuxiliaryFiles => false,
	 Authors => {{Name => "David Eisenbud", Email => "de@msri.org"},
	    {Name => "Michael Perlman", Email => "mperlman@umn.edu"}, 
	    {Name => "Ritvik Ramkumar", Email => "ritvikr@cornell.edu"},
	    {Name => "Deepak Sireeshan"},
	    {Name => "Aleksandra Sobieska", Email => "asobieska@wisc.edu"},
	    {Name => "Teresa Yu", Email => "twyu@umich.edu"},
	    {Name => "Jacob Zoromski", Email => "jzoromsk@nd.edu"} },
	 Headline => "Creating random ideals for experiments",
	 PackageExports => {"DGAlgebras", "RandomIdeals", "MonomialOrbits"},
         DebuggingMode => false )
export {
    "mixedIdeal",
    "mixedIdealList",
    "randomMonomialIdealHilbertFunction"
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

randomMonomialIdealHilbertFunction = method()
randomMonomialIdealHilbertFunction (Ring, List) := MonomialIdeal => (S, hf) -> (
    idealsL := hilbertRepresentatives(S, hf);
    ranN := random(0,#hf-1);
    I := idealsL#ranN;
    I + ((ideal vars  S)^(#hf + 1))
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

doc ///
Key
  randomMonomialIdealHilbertFunction
Headline
 Gives a random Artin monomial ideal with a given Hilbert function.
Usage
  I  = randomMonomialIdealHilbertFunction(S,hf)
Inputs
 S: Ring
 hf: List
Outputs
 I: MonomialIdeal
Description
  Text
    This function provides a random monomial ideal given a specified Hilbert function $\text{HF}$  in the
    form of a list $\{\text{HF}(1),\ldots,\text{HF}(d)\}$. It is assumed that $\text{HF}(n)=0$ for
    $n>d$.
  Example
    S = ZZ/101[a,b,c,d]
    hf = {4,7,10}
    I = randomMonomialIdealHilbertFunction(S, hf)
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


