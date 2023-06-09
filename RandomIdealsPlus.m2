newPackage(
         "RandomIdealsPlus",
         Version => "0.1",
         Date => "",
	 AuxiliaryFiles => false,
	 Authors => {{Name => "David Eisenbud", Email => "de@msri.org"},
	    {Name => "Michael Perlman", Email => "mperlman@umn.edu"}, 
	    {Name => "Ritvik Ramkumar", Email => "ritvikr@cornell.edu"},
	    {Name => "Deepak Sireeshan", Email => "dsbx7@umsystem.edu"},
	    {Name => "Aleksandra Sobieska", Email => "asobieska@wisc.edu"},
	    {Name => "Teresa Yu", Email => "twyu@umich.edu"},
	    {Name => "Jacob Zoromski", Email => "jzoromsk@nd.edu"} },
	 Headline => "Creating random ideals for experiments",
	 PackageExports => {"DGAlgebras", "RandomIdeals", "MonomialOrbits", "Macaulay2Doc"},
         DebuggingMode => false )
export {
    "mixedIdeal",
    "mixedIdealList",
    "randomMonomialIdealHilbertFunction",
    "randomNnomialIdeal",
    "randomMixedIdeal",
    "Pure",
    "randomArtinMonomialIdeal", 
    "Randm"
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

randomNnomialIdeal = method(Options => {Pure => false})
randomNnomialIdeal(ZZ,List,Ring) := o -> (N,Degs,S) -> (
    if o.Pure == false then coeffs := ultimate (coefficientRing, S) else coeffs = {1,-1};
    ideal apply(Degs, d -> (
	    isNew := false; --whether the generated monomial is in the support of the new generator already
	    newGen := 0_S;
	    numTerms := 0;
	    while numTerms < N do (
	        while isNew  == false do (
		    m := randomMonomial(d,S);
		    isNew = all(flatten entries monomials(newGen), i -> i != m);
		    );
		if o.Pure == false then newGen = newGen + random(coeffs)*m else newGen = newGen + (random(coeffs))_0*m;
		isNew = false;
		numTerms = numTerms + 1;
		);
	    newGen
	    )
	)
    )

randomMixedIdeal = method(Options => {Pure => false})
randomMixedIdeal(List,List,Ring) := o -> (TermList,DegList,S) -> (
    sum apply(#DegList, i -> randomNnomialIdeal(TermList#i,DegList#i,S,Pure => o.Pure))
    )

randomArtinMonomialIdeal = method(Options => {Randm => false})
randomArtinMonomialIdeal(List, PolynomialRing) := Ideal => o -> (Dlis, A) -> (
     if o.Randm == true then (
     	 I := trim(randomMonomialIdeal(Dlis,A));
     	 indet := gens A;
     	 K := first entries gens I;
     	 m :=(#indet) - 1;
     	 if  dim(A/I) != 0 then (
	     leastcm := lcm(K);
	     L :={};
	     for i from 0 to m do (
	     	 deg := (degree(indet_i,leastcm)) + 1;
	     	 L = append(L,(indet_i)^deg);
	     	 );
	     newI := trim(ideal(join(K,L)));
	     return newI
	     )
     	 else return I
	 )  
     else if (#Dlis != #gens A) then error "Use option 'rand' or choose appropriate list"
     else if (#Dlis == #gens A) then  (
	 indet = gens A;
    	 n := (#indet)-1;
    	 mir := max(Dlis);
    	 mar := sum(Dlis) - n + 1;
    	 detList := for i from 0 to n list (indet_i)^(Dlis_i);
    	 total := random(2,5);
    	 for i from 1 to total do (
	     criteria1 := false;
	     while criteria1 == false do (
	    	 mdeg := random(mir,mar);
	    	 mmon := randomMonomial(mdeg, A);
	    	 criteria1 = (mmon % ideal(detList)) != 0_A;
	    	 );
	     detList = append(detList,mmon);
    	     );
    	 return trim(ideal(detList)) 
	 ) 
)	
     
randomArtinMonomialIdeal (List, ZZ, PolynomialRing) := Ideal => o -> (Dlis, d, A) -> (
    indet := gens A;
    n := (#indet)-1;
    mir := max(Dlis);
    mar := sum(Dlis) - n + 1;
    detList := for i from 0 to n list (indet_i)^(Dlis_i);
    for i from 1 to d do (
	criteria2 := false;
	while criteria2 == false do (
	    mdeg := random(mir,mar);
	    mmon := randomMonomial(mdeg, A);
	    criteria2 = (mmon % ideal(detList)) != 0_A;
	    );
	detList = append(detList,mmon);
    );
    return trim(ideal(detList)) 
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
  Collection of functions  to create random ideals (or) a list of random ideals of Monomial, mixed or Artin Ideals.
///


doc ///
Key
  mixedIdeal
  (mixedIdeal, Ring, List, List) 
Headline
 Randomly generates an ideal which is a sum of a monomial and a binomial ideal. 
Usage
  I = mixedIdeal(R,L1,L2)   
Inputs
 R: Ring
  Polynomial Ring
 L1: List
  Degree of the binomials of $R$  that are part of the generating set of the ideals.
 L2: List
  Degree of the monomials of $R$ that are part of the generating set of the ideals.
Outputs
 I: Ideal
Description
  Text
    Mixed ideal(in this context) is the sum of a monomial ideal and a binomial ideal. The function
    produces a random mixed ideal. 
  Example
    R = ZZ/101[a,b,c]
    Blist = {3,4}
    Mlist = {2,3}
    I = mixedIdeal(R, Blist, Mlist)
///


doc ///
Key
  mixedIdealList
  (mixedIdealList, Ring, ZZ, List, List)
Headline
 Randomly generates a list of ideals which are sums of monomial and binomials ideals. 
Usage
  L  = mixedIdealList(R, i, L1, L2)   
Inputs
 R: Ring
  Polynomial ring.
 i: ZZ
  Number of ideals you want. 
 L1: List
  Degree of the binomials of $R$  that are part of the generating set of the ideals. 
 L2: List
  Degree of the monomials of $R$ that are part of the generating set of the ideals.  
Outputs
 L: List
  List of ideals that are generated by random monomials and binomials. 
Description
  Text
    Mixed ideal(in this context) is the sum of a monomial ideal and a binomial ideal. The function
    produces a list of random mixed ideals. 
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



doc ///
Key
  randomNnomialIdeal
  (randomNnomialIdeal, ZZ, List, Ring)
  [randomNnomialIdeal, Pure]
Headline
 Computes a random ideal generated by homogeneous elements with a fixed number of terms.
Usage
  I  = randomNnomialIdeal(N,Degs,S)   
Inputs
 N: ZZ
 Degs: List
 S: Ring
Outputs
 I: Ideal
Description
  Text
    This function produces a random monomial, binomial, trinomial, etc., ideal with given degrees.
    For example, we can produce a random trinomial (N = 3) ideal with two generators of degrees 2 and 3
  Example
    S = ZZ/101[a,b,c]
    N = 3
    Degs = {2,3}
    I = randomNnomialIdeal(N,Degs,S)
  Text
    With the option Pure => true, the coeffecients in each generator will be 1 or -1.
  Example
    S = ZZ/101[a,b,c]
    N = 4
    Degs = {4}
    I = randomNnomialIdeal(N,Degs,S, Pure => true)
///


doc ///
Key
  randomMixedIdeal
  (randomMixedIdeal, List, List, Ring)
Headline
 Computes a random ideal generated by homogeneous elements with specified numbers of terms and specified degrees.
Usage
  I  = randomMixedIdeal(TermList,DegList,S)   
Inputs
 TermList: List
 DegList: List
 S: Ring
Outputs
 I: Ideal
Description
  Text
    This function gives a random ideal whose generators have fixed numbers of terms specified by TermList.
    For example TermList = {2,3} means the ideal will be generated by binomials and trinomials.
    DegList is a list of lists, the degrees of the generators of with the respective number of terms.
    For example, to generate an ideal with one binomial of degree 3 and two trinomials of degrees 3 and 4, 
    we do the following:
  Example
    S = ZZ/101[a,b,c]
    TermList = {2,3}
    DegList = {{3},{3,4}}
    I = randomMixedIdeal(TermList,DegList,S)
  Text
    The option Pure => true makes all of the coefficients either 1 or -1.
  Example
    S = ZZ/101[a,b,c]
    TermList = {2,3}
    DegList = {{5},{3}}
    I = randomMixedIdeal(TermList,DegList,S,Pure => true)
///

doc ///
Key
 randomArtinMonomialIdeal
 (randomArtinMonomialIdeal, List, PolynomialRing)
 (randomArtinMonomialIdeal, List, ZZ, PolynomialRing)
Headline
 Randomly generates an Artin Monomial Ideal. 
Usage
 I = randomArtinMonomialIdeal(ranList, R)
 J = randomArtinMonomialIdeal(detList, n, R)
 K = randomArtinMonomialIdeal(detList, R) 
Inputs
 ranList: List 
  List of degrees of monomials of $R$ that are part of the generating set of the ideal. 
 R: PolynomialRing
 detList: List
  List of powers of generators of $R$ that you want in the ideal generated.  
 n: ZZ
  The number of monomials adjoined to the generating set of ideals such that each monomial added
  to the list is not in the ideal generated by the preceding list. 
Outputs
 I: Ideal
  Artin Ideal generated by random list of monomials. 
 J: Ideal
  Artin Ideal generated by the powers of generators of the maximal ideal along with $n$ other
  random  monomials.
 K: Ideal
  Artin Ideal generated by the powers of generators of the maximal ideal along with a
  random number, $2 \leq n \leq 20$  of random monomials.  
Description
  Text
   The function can be used in one of these three ways. 
   The following example computes a random Artin Monomial ideal where the user needs to 
   provide the degrees of monomials they need in its set of generators. 
  Example
   R = ZZ/101[a,b,c]
   DegLis = {2,3}
   I = randomArtinMonomialIdeal(DegLis, R, Randm => true)
  Text
   In the example above, rand is an optional imput which is false by default. The two monomials
   of any degree cannot generate an Artin ideal. The function takes these monomials and adds
    appropriate powers of generators of maximal ideal $(a,b,c)$ such that 
    it is the smallest Artin monomial ideal that contains the random monomials 
    the function picked. 
  Example
   R = ZZ/101[a,b,c]
   DegLis = {4,7,9}
   I = randomArtinMonomialIdeal(DegLis, R)
  Text
   In this example, we use the default value of the options. The function randomly 
   picks the number of extra irredundant monomials to be added provided and considers Dlist as
   the degrees of the interminates of the ring. S0, {4,7,9} indicates that the monomials 
   {$a^4, b^7, c^9$} are added to the generating list. It is important to notice that if the
   length of the input list is smaller than the number of indeterminates, then the function
   displays an error message. 
  Example
   R = ZZ/101[a,b,c]
   DegLis = {4,7,9} 
   I = randomArtinMonomialIdeal(DegLis, 3, R)
  Text
   In this last example, the input list {4,7,9} indicates the power of indeterminates as before.
   The generating set will also have 3 irredundant monomials added to it such that 
   if $m_i$ is the $i$th addition, then $m_i \notin (a^4, b^7, c^9, m_0, \cdots, m_{i-1})$.
  
///  

-* Test section *-
TEST /// -* [insert short title for this test] *-
setRandomSeed 0;
S = ZZ/101[a..d];
hf = {4,7,10};
I = randomMonomialIdealHilbertFunction(S,hf)
assert(I == ideal "a2, ab, b2, a4, a3b, a3c, a3d, a2b2, a2bc, a2bd, a2c2, a2cd, a2d2, ab3, ab2c, ab2d, abc2, abcd, abd2, ac3, ac2d, acd2, ad3, b4, b3c, b3d, b2c2, b2cd, b2d2, bc3, bc2d, bcd2, bd3, c4, c3d, c2d2, cd3, d4")
setRandomSeed 0
assert(randomNnomialIdeal(3,{3},S) == ideal( -30*a^2*b + 24*b^2*c + 19*a*c*d))
setRandomSeed 0
assert(randomMixedIdeal({2,3},{{2},{2}},S) == ideal(-10*a*d - 8*c*d, 19*a*c + 39*b*c - 24*c*d))
R = ZZ/101[a,b,c]
D = {2,3,4}
E = {2,3}
setRandomSeed(0)
assert(randomArtinMonomialIdeal(D,R) == ideal(a^2, b^3, c^4, b*c^3, a*c^3, b^2*c^2, a*b*c^2, a*b^2*c))
setRandomSeed(0)
assert(randomArtinMonomialIdeal(E, R, Randm => true) == ideal(a,c^2,b^3))
setRandomSeed(0)
assert(randomArtinMonomialIdeal(D,2,R) == ideal(a^2, b^3, c^4, b*c^3, b^2*c^2)
setRandomSeed(0) 
assert(mixedIdeal(R,D,E) == ideal(a^2-b^2, a*b^2-b^3, b^3*c - a*b*c^2, c^2, b^2*c))
setRandomSeed(0)   
assert(mixedIdealList(R,1,D,E) == {ideal(b*c, a*b +21*b^2, a*c^2, b^3)})    
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


