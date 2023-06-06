newPackage(
	"SocleSummands",
	Version => "0.3",
	Date => "June, 2023",
	AuxiliaryFiles => false,
	Authors => {{Name => "David Eisenbud", Email => "de@msri.org"},
	    {Name => "Hai Long Dao"}},
	Headline => "Experiments on infinite resolutions",
	PackageExports => {"DGAlgebras", "MonomialOrbits","RandomMonomialIdeals"},
	DebuggingMode => true)
export {
    "summandExamples",
    "socle",
    "socleSummands",
    "socleSummandsSemigroup",
    "isBurch",
    "burchIndex",
    --
    "Boundaries"

}
    


 summandExamples = (S = ZZ/101[a,b,c];
    hashTable{(1, ideal"a4,b4,c4,ab3"),
    (2,ideal"a4,b4,c4,ab3,cb3"),
    (3, ideal"a4,b4,c4,ab3,bc3"),
    (4, ideal"a4,b4,c4,abc"),
    (5, ideal"a4,b4,c4,ab3, b2c2"),
    (6, ideal"a4,a3b,b3,abc2,c5")}
)

socle = method()
socle Ring := R -> (
    (ideal 0_R):(ideal gens R)
    )
socle Module := M -> (
    R := ring M;
    (0_R*M):(ideal gens ring M)
    )


socleSummands = method(Options => {Verbose => false, Boundaries => true})

-*
socleSummands Module := o-> M -> (
   mm := ideal gens ring M;
   ss := socle M;
    if o.Verbose == false then
       numcols compress (gens ss % (mm*M)) 
          else
       flatten degrees source (gens ss % (mm*M))
    )
*-

socleSummands ChainComplex := o -> C -> (
    --for i from 1+min C to max C, return
    --numgens  (socle of target C.dd_i)/(m*image C.dd_i)
    --or the degrees of these gens (slower)
    R := ring C;
    mm := ideal gens R;
    socR := trim (ideal(0_R):mm);
    if o.Boundaries == true then
       for i from min C + 1 to max C list 
           imageSocleSummands(C.dd_i, socR, Verbose => o.Verbose)
    else for i from min C + 1 to max C list (
	   if i == min C + 1 then imageSocleSummands(C.dd_i, socR,  Verbose => o.Verbose)
	   else imageSocleSummands(syz C.dd_(i-1), socR, Verbose => o.Verbose)
	   ))

socleSummands(Module,ZZ) := o -> (M,ell) -> (
    F := res(M, LengthLimit => ell);
    socleSummands (F,o)
    )

imageSocleSummands = method(Options =>{Verbose => false})
imageSocleSummands(Matrix,Ideal) := o -> (phi, socR) ->(
    --phi:F \to G map of free modules
    --socR = gens trim ideal(0_R):mm;
    --the ideal of the socle in ring F
    --compute the number and optionally degrees of 
    --the (socle of target phi)/(m*im phi)
    
    --ss is a non-minimal map to source G**socR, a free R-module
    --whose basis corresponds to the socle of G. Thus numrows ss = rank socle G.
    --the map ss surjects to the summand
    --corresponding to socle elements in mm*(image phi), 
    --thus rank (Rbar**ss) is the rank of the summand of the socle of G,
    --contained in mm(image phi).
    --the cokernel of Rbar**ss is thus the vector space of
    --socle summands.
    
    R := ring phi;
    F := source phi;
    G := target phi;
    Rbar := R/(ideal vars R);
    ss := modulo(G**gens socR, phi**vars R);
      socSummands := coker (Rbar**ss);
     if not o.Verbose then degree socSummands else
          degrees socSummands
    )

-*
socleSummands = method(Options => {Verbose => false})
socleSummands(ChainComplex, ZZ) := o -> (C,m) -> (
    --for i from 1+min C to m, return
    --numgens  (socle of target C.dd_i)/(m*image C.dd_i)
    --or the degrees of these gens (slower)
    R := ring C;
    mm := ideal gens R;
    socR := gens(ideal(0_R):mm);
    Rbar := R/mm;
    for i from min C + 1 to m list(
	phi := C.dd_i;
        ss := modulo((target phi)**socR, phi**vars R);
	    if o.Verbose == false then
        numrows ss - rank(Rbar**ss)
          else
       (degrees gens prune coker ss)_0
       )
)
*-

hasSocleSummand = method()
hasSocleSummand Module := M -> (
    0 != socleSummands M)


inSemigroup = method()
inSemigroup(ZZ,List) := (a,L) -> (
	val := false;
	if (#L == 1) then (
		val = (a%(L#0) == 0)
		);
	if (#L >= 2) then (
		for i from 0 to a//(L#0) do (
			val = val or inSemigroup(a-i*(L#0),L_{1..(#L-1)});
			if val then break;
			);
		);
	val
	)

socleSummandsSemigroup = method(Options => {Verbose => false, Boundaries => true})
socleSummandsSemigroup ChainComplex :=  o -> F -> (
    L := socleSummands(F,o);
    P := positions(L, ell -> ell != 0);
    	gens := {0}; --identity element of the semigroup
	for p in P do
	        if  (not inSemigroup(p+1,gens)) then (gens = gens | {p+1});
	gens
	)

socleSummandsSemigroup(Ideal, ZZ) := o -> (I,n) -> (
    R := ring I/I;
    F := res(coker vars R, LengthLimit => n);
    socleSummandsSemigroup(F,o)
    )


isBurch = method()
isBurch Ideal := I -> (
    R := ring I;
    mm := ideal gens R;
    mm != (mm*I):(I:mm)
    )

burchIndex = method()
burchIndex Ideal := I -> (
    R := ring I;
    mm := ideal gens R;
    degree(mm/((mm*I):(I:mm)))
    )


      -* Documentation section *-
      
beginDocumentation()


doc ///
Key
 "summandExamples"
Headline
 Creates hashtable of monomial ideal examples
Usage
 H = summandExamples ()
Outputs
 H: HashTable
Description
  Text
    Produces a hashTable of monomial ideals in ZZ/101[a,b,c]
  Example
    summandExamples
///

doc ///
Key
 socle
 (socle, Ring)
 (socle, Module)
Headline
 Computes the socle of a ring or module
Usage
 I = socle R
 N = socle M
Inputs
 R: Ring
 M: Module
Outputs
 I: Ideal
 N: Module
Description
  Text
   The socle of a ring or module is
   the sum of the minimal nonzero ideals or submodules
  Example
   S = ZZ/101[a,b]
   I = (ideal(gens S))^2
   R = S/I
   socle R
   M = coker random(R^2, R^{-1})
   socle M
   degree socle M
///

doc ///
Key
 burchIndex
 (burchIndex, Ideal)
Headline
 Computes the Burch index of an ideal
Usage
 bi = burchIndex(I)
Inputs
 I: Ideal
Outputs
 bi: ZZ
Description
  Text
   If $(S,n,k)$ is a regular local ring with
   ideal $I$ such that $R=S/I$ is a depth 0
   local ring with residue field $k$, then the
   Burch index of $I$ is the dimension of the
   vector space $n/(In:(I:m))$.
   The ring $R$ is called Burch if the Burch index
   of $I$ is at least $1$. If the Burch index of $I$
   is at least $2$, then the minimal free resolution of a non-free
   $R$-module $M$ has the residue field as a direct
   summand of the $n$th syzygy module of $M$ for $n\geq 7$.
  Example
   S = ZZ/101[a,b]
   I = (ideal(gens S))^2
   burchIndex(I)
///

doc ///
Key
 isBurch
 (isBurch, Ideal)
Headline
 Determines whether an ideal is a Burch ideal
Usage
 isBur = isBurch(I)
Inputs
 I: Ideal
Outputs
 isBur: Boolean
Description
  Text
   An ideal $I$ of a local ring $(S,m,k)$ is an ideal where $(m*I):(I:m)$
   is strictly contained in $m$. Equivalently, $I$ is Burch if and only if
   the second syzygy of $k$ over $S/I$ contains $k$ as a direct summand.
  Example
   S = ZZ/101[a,b]
   I = (ideal(gens S))^2
   isBurch(I)
   J = ideal(a^2,b^2)
   isBurch(J)
///





doc ///
Key
 socleSummands
 (socleSummands, Module, ZZ)
 (socleSummands, ChainComplex)
 [socleSummands, Verbose]
 [socleSummands, Boundaries]
Headline
 Computes the number of socle summands of the images of the differentials of a chain complex or free resolution of a module
Usage
 S = socleSummands (M,n)
 S = socleSummands C
Inputs
 M: Module
 n: ZZ
 C: ChainComplex
Outputs
 S: List
Description
  Text
   Given a chain complex, this function calculates the number of socle summands in the images of each differential.
   If one instead inputs a module M and an integer n, the function will calculate the socle summands in the images of the differentials of the free resolution of M truncated after n steps. 
   Note that function ignores the 0-th differential of the resolution.
   One may use the option Boundaries => false to instead find the number of socle summands of the cycles of the complex (which will give the same output for a free resolution)
   One may use the option Verbose => true to instead find the degrees of the socle summands
  Example
   S = ZZ/101[a,b,c]
   I = ideal(a^4,b^4,c^4,a*b^3,b*c^3)
   R = S/I
   mm = ideal gens R
   F = res( mm, LengthLimit => 6)
   socleSummands F
  Text
    We will get the same ouput if we instead input the module and length limit:
  Example
    S= ZZ/101[a,b,c]
    I = ideal(a^4,b^4,c^4,a*b^3,b*c^3)
    R = S/I
    mm = ideal gens R
    socleSummands((R^1/mm),6)
  Text
    We can instead use the cycles of the free resolution (which again will give the same output):
  Example
    S= ZZ/101[a,b,c]
    I = ideal(a^4,b^4,c^4,a*b^3,b*c^3)
    R = S/I
    mm = ideal gens R
    socleSummands((R^1/mm),6, Boundaries => false)
///


doc ///
Key
 socleSummandsSemigroup
 (socleSummandsSemigroup, ChainComplex)
 (socleSummandsSemigroup, Ideal, ZZ)
Headline
 Computes the semigroup generated by the socle summands of the boundaries or cylces of a chainComplex
Usage
 L = socleSummandsSemigroup F
 L = socleSummandsSemigroup(F, Boundaries => False)
 L = socleSummandsSemigroup(I,n)
Inputs
 F: ChainComplex
 I: Ideal
 n: ZZ
Outputs
 L: List
Description
  Text
   The homological degrees of the socle summands of the 
   boundaries or cycles of a chainComplex form a semigroup. 
   If the input is a chain complex F, socleSummandsSemigroup 
   computes the socle summands of the boundaries of F. 
   To compute the socle summands of the cycles
   set Boundaries => False. 
  
   socleSummandsSemigroup of an ideal computes the semigroup 
   of the residue field. 
  Example
   S = ZZ/101[a,b,c]
   I = ideal"a4,a3b,b3,abc2,c5"
   socleSummandsSemigroup(I,7)
  Example
   S = ZZ/101[a,b,c]
   I = ideal "a4,a3b,b3,a2c2,c5"
   socleSummandsSemigroup(I,7)
  Example
   S = ZZ/101[a,b,c]
   I = ideal "a4,a3b,b3,a2c2,c5"
   R = S/I
   M = cokernel random(R^1,R^{-1}) 
   F
   socleSummandsSemigroup(F)
  Example
   S = ZZ/101[a,b,c] 
   I = ideal "a4,a3b,b3,a2c2,c5"
   R = S/I
   F = koszul(matrix{{a,b,c}})
   socleSummandsSemigroup(F, Boundaries => false)

///


-* Test section *-
TEST///
L = for i from 1 to 4 list(
I = summandExamples#i;
socleSummandsSemigroup(I,7)
)
assert (L == {{0}, {0,2,3}, {0,3,4,5}, {0,4,5,6,7}})
///


end--

uninstallPackage "SocleSummands"
restart
installPackage "SocleSummands"


restart
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"
needsPackage "RandomIdeals"

setRandomSeed 0



for i from 1 to 4 do(
<<(I = summandExamples#i)<<" "<<socleSummandsSemigroup(I, 7)<<endl;
)
koszul vars ring I
socleSummands oo
I = summandExamples#1
socleSummandsSemigroup(I,7)
F = res(coker vars (ring I/I), LengthLimit => 7)
socleSummandsSemigroup F
socleSummands F
S = ZZ/101[a..c]
I = mixed (S,{3,4,4,5},{2,3, 4,5,6, 7})
dim I
R = S/I
F = res (coker vars R, LengthLimit => 7)
socleSummandsSemigroup F
socleSummandsSemigroup(I, 5)
isBurch I
burchIndex I



restart
load "deExplore.m2"
S = ZZ/101[a..c]
L = for i from 1 to 1000 list (
    I := trim mixed(S, {3,4},{4,5,6,7});
    if dim I> 0 then continue;
    I
	);
#L
Lnb = select(L, I -> not isBurch I);#Lnb
Lng = select(Lnb, I -> not isGolod (S/I));#Lng
tally apply(Lng, I-> numgens I)
elapsedTime for I in Lng list kSS(I, 7)

R = S/Lng_4
F = res (coker vars R, LengthLimit => 7)
socleSummands F




use S
S === ring Lng_3
R = (ring Lng_3)/(Lng_3)
F = res(coker vars R, LengthLimit => 7)

socleSummands F
for i from 1 to length F do(
     <<socleSummands(prune image F.dd_i);
	 <<endl;
     )
elapsedTime ll = apply(Lng, I -> (kSS(I, 7)))
gens prune((socle image F.dd_3))%((ideal vars R)*image F.dd_3)
(socle R^1)*(image F.dd_3)
ideal socle (R^1)
mm = ideal vars R
socR = (ideal 0_R):mm
socF = socR*F_3
(gens socF) %(mm*image F.dd_3)

R = S/Lng_1
betti (F = res (coker vars R, LengthLimit => 7))
use S
I = Lng_1
isBurch I
mm = ideal vars S
BI = (mm*I):(I:mm)
prune image F.dd_2

tally(Lng/dim)
LL = for I in Lng list (I + ideal random(S^1, S^{-1}))
LL = select(LL, I-> dim I == 0 and not isBurch I and not isGolod (S/I))

kSS(LL_0, 10)
use (R = S/LL_0)
betti res(coker vars R, LengthLimit => 10)

uninstallPackage "MonomialOrbits"
installPackage "MonomialOrbits"

peek loadedFiles
p
debug needsPackage "MonomialOrbits"
code methods orbitRepresentatives
S = ZZ/101[a..c]
code methods toLis
orbitRepresentatives(S,{0,3})

