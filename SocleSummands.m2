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
    "Examples",
    "socle",
    "socleSummands",
    "socleSummandsSemigroup",
    "isBurch",
    "burchIndex",
    --

}
    


 Examples = (S = ZZ/101[a,b,c];
    hashTable{(1, ideal"a4,b4,c4,ab3"),
    (2,ideal"a4,b4,c4,ab3,cb3"),
    (3, ideal"a4,b4,c4,ab3,bc3"),
    (4, ideal"a4,b4,c4,abc"),
    (5, ideal"a4,b4,c4,ab3, b2c2")}
)

socle = method()
socle Ring := R -> (
    (ideal 0_R):(ideal gens R)
    )
socle Module := M -> (
    R := ring M;
    (0_R*M):(ideal gens ring M)
    )

socleSummands = method(Options => {Verbose => false})
socleSummands Module := o-> M -> (
    mm := ideal gens ring M;
    ss := socle M;
    if o.Verbose == false then
       numcols compress (gens ss % (mm*M)) 
          else
       flatten degrees source (gens ss % (mm*M))
    )

socleSummands ChainComplex := o -> C -> (
    --for i from 1+min C to max C, return
    --numgens  (socle of target C.dd_i)/(m*image C.dd_i)
    --or the degrees of these gens (slower)
    R := ring C;
    mm := ideal gens R;
    socR := trim (ideal(0_R):mm);
    for i from min C + 1 to max C list 
        imageSocleSummands(C.dd_i, socR, o)
       )

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

socleSummandsSemigroup = method()
socleSummandsSemigroup ChainComplex :=  F -> (
    L := socleSummands F;
    P := positions(L, ell -> ell != 0);
    	gens := {0}; --identity element of the semigroup
	for p in P do
	        if  (not inSemigroup(p+1,gens)) then (gens = gens | {p+1});
	gens
	)

socleSummandsSemigroup(Ideal, ZZ) := (I,n) -> (
    R := ring I/I;
    F := res(coker vars R, LengthLimit => n);
    socleSummandsSemigroup F
    )

isBurch = I -> (
    R := ring I;
    mm := ideal gens R;
    mm != (mm*I):(I:mm)
    )
burchIndex = I -> (
    R := ring I;
    mm := ideal gens R;
    degree(mm/((mm*I):(I:mm)))
    )


      -* Documentation section *-
      
beginDocumentation()

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



-* Test section *-
TEST///
L = for i from 1 to 4 list(
I = Examples#i;
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

mixedIdeal = (S, BList, MList) -> (randomPureBinomialIdeal(BList, S)+randomMonomialIdeal(MList,S))   

for i from 1 to 4 do(
<<(I = Examples#i)<<" "<<socleSummandsSemigroup(I, 7)<<endl;
)
koszul vars ring I
socleSummands oo
I = Examples#1
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

