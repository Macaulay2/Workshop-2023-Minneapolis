newPackage(
	"SocleSummands",
	Version => "0.2",
	Date => "January, 2022",
	AuxiliaryFiles => false,
	Authors => {{Name => "David Eisenbud", Email => "de@msri.org"},
	    {Name => "Hai Long Dao"}},
	Headline => "Determining (kSi) Property in Artin Local Rings",
	PackageExports => {"Points", "DGAlgebras", "MonomialOrbits"},
	DebuggingMode => true)

export {
        "Examples",
	"socle", -- compute socle(M) as the image of a map to M
	"hasSocleSummand", -- this version works even when the module is not a syzygy.
	"hasSocleSummand1", -- this version works only when the module is given as a sub of a free with the same socle.
	"findSocleSummands1", -- this version works only when the module is given as a sub of a free with the same socle.	
	"socleSummands", -- now applies to modules and resolutions and ideals. Takes the place of kSI etc
    	"weakLin",
	"ent1",--entries of the maps in a chain cplx
        "minimalSummand",
--	"kSi", -- find the terms in the free res of k with socle summands
--	"kSRes",-- same, assuming the complex it gets is acyclic
--	"kkSi",
	"kSS", -- outputs the semigroup generators, where the socle summands in res vars R are.
--	"kSS2",	
--	"kSI",
--	"MSi",
--	"MSiList",
	"linearFormsInSyzygy",
--	"lin1",
--    	"koSi",
        "socleSummandGeneration", -- version for a DGA needs work
	"koS",
	"TorMap",
	"TorMaps",
	"isTorSurjective",
	"nextList",
	"ListToMultiDegree",
	"monomialList",
	"monomialIdealList",
	"isReducedList",
	"test1",
	"golodTests",
	"TorMapk",
	"TorMapConjecture",
	"lin7",
	"dist",
	"bm",
	"burchIndex",
	"isBurch",
	"linearFormsInRes",
	"KtoR", -- constructs first stages of the Tate res of k as a DGAlgebra.
--
--
    	"ringOfRandomPoints",
	"resOfRandomPoints",
	"randomMonomial",
	"randomArtinIdeal",
	"randomArtinDegreeIdeal",
	"randomArtinMonomialIdeal",
	"randomArtinDegreeMonomialIdeal",
	"randomArtinRing",
	"randomArtinDegreeRing",
	"randomArtinMonomialRing",
	"randomArtinDegreeMonomialRing",
	"ran",
	"Check", -- option for ran
	"Bound",
	"Numinfo", -- option for hasSocleSummand1
	"Nonminimal",
	"Resolution" -- option for KtoR
	}
///
--todo: remove documention of archaic functions
restart
uninstallPackage "SocleSummands"
restart
installPackage "SocleSummands"
loadPackage "SocleSummands"
viewHelp SocleSummands
///

koSi = method()
koSi Ring := R -> (
    K := koszul vars R;
    mm := ideal vars R;
    socleR := ideal(0_R):mm;
    apply (length K, i -> (
	    N := socleR*K_(i+1); 
--    	    (gens N) % (mm*image syz K.dd_(i+1))
	    (gens N) % map(target gens N,,
		          (vars R)**syz K.dd_(i+1))
    	))
)


socleSummandGeneration = method()
-- test whether the socle summands of homological degree j+2
-- in the Koszul complex are generated,
-- as a module over the exterior algebra,
-- by the ones in degree j+1, for j = 0 to numvars -2.
socleSummandGeneration Ring := List => R -> (
n := numgens R;
mmR := ideal vars R;
socR := (ideal 0_R):mmR;
K := koszulComplexDGA R;
L := toComplex K;
A := K.natural;
soc := apply(n-1, i-> gens (socR*L_(i+1)));
cyc := apply(n-1, i-> syz L.dd_(i+1));
socSummands := apply(n-1, i-> soc_i % map(target soc_i,,
	                           (vars R**cyc_i)));
s := apply (n-2, j ->(
sum apply(n,i-> image(
	             (dgAlgebraMultMap(K,A_i))_(j+1)*socSummands_j
	     ))
     ));
apply(n-2, j->(
(socSummands_(j+1) % gens s_(j) == 0, 
 numgens prune coker socSummands_(j+1))))
)

///
restart
loadPackage("SocleSummands", Reload => true)
loadPackage("DGAlgebras", Reload => true)
R = ZZ/101[a,b,c]/ideal"a3,abc,b2,c2"
isBurch ideal R
socleSummandGeneration R
K = koszulComplexDGA R
socleSummandGeneration(K,3)
KK = KtoR(R,2)

socleSummandGeneration(KK_2,3)
i=0;j=0
A
socSummands
             (dgAlgebraMultMap(K,A_i))_(j+1)*socSummands_j
K = KK_2
///

socleSummandGeneration(DGAlgebra,ZZ) := List => (K,n) -> (
R := K.ring;

mmR := ideal vars R;
socR := (ideal 0_R):mmR;
L := toComplex(K,n);
A := K.natural;
soc := apply(n-1, i-> gens (socR*L_(i+1)));
cyc := apply(n-1, i-> syz L.dd_(i+1));
socSummands := apply(n-1, i-> soc_i % (mmR*(image L.dd_(i+2))));
s := apply (n-2, j ->(
	sum apply(numgens A, i-> image(
	             (dgAlgebraMultMap(K,A_i))_(j+1)*socSummands_j
	     ))
     ));
apply(n-2, j->(
(socSummands_(j+1) % gens s_(j) == 0, 
numgens prune coker socSummands_(j+1))))
)


///
needsPackage "SocleSummands"
S = ZZ/101[a_1..a_5]
mS = ideal gens S
R = S/(mS^3)
socleSummandGeneration R
///
minimalSummand = method()
minimalSummand Matrix := Matrix => M ->(
    --decompose a homogeneous M into a trivial map and a minimal summand;
    --return the minimal summand
    Mm := coker M;
    presentation prune Mm
    )

ent1 = method()
ent1 ChainComplex := F ->(
    ell := length F;
    netList for i from 1 to ell list trim ideal F.dd_i
    )

ent1 Matrix := Ideal => M -> (
    trim ideal minimalSummand M
    )
///
restart
loadPackage "SocleSummands"
S = ZZ/101[a,b]
M = map(S^1, S^{2:-1}, matrix"a,b")++map(S^1, S^1, 1)
M' = random(target M, target M) *M* random(source M, source M)
ent1 M'
///

bm = method()
bm Ideal := Ideal => I -> (
    S := ring I;
    mm := ideal vars S;
    mm*I:(I:mm)
    )
burchIndex = method()
burchIndex Ideal := ZZ => I ->(
    mm := ideal vars ring I;
    numgens prune (mm/bm(I))
    )    
isBurch = method()
isBurch Ideal := Boolean => I -> burchIndex I > 0

linearFormsInSyzygy = method()
linearFormsInSyzygy(Module, ZZ) := (N,n) -> (
    S := ring N;
    mm2 := (ideal vars S)^2;
    F := res(N, LengthLimit => n);
    apply(1+length F, i -> 
	       numcols compress((gens trim ideal F.dd_i)%mm2)
	       )
	   )
linearFormsInSyzygy(Ideal, ZZ) := (N,n) -> linearFormsInSyzygy(module N,n)
///
S = ZZ/101[a,b,c]
I = ideal"a3,b3,c3"*ideal(a,b,c)
R= S/I
J = module ideal"abc"
linearFormsInSyzygy(J, 4)
///

--(i, true) means that the entries of F.dd_(i+1) are contained in the sum of the previous.
weakLin = method()
weakLin (Module, ZZ) := (M, n) ->(
    S = ring M;
    F := res(M, LengthLimit => n);
    I := trim ideal F.dd_1;
    Jtot := I;    
    for i from 1 to n do(
	I = trim ideal F.dd_(i+1);
	t:= (gens I % Jtot == 0);
	<<(i,t)<<endl;
    	Jtot = trim(Jtot+I))
    )

weakLin(Ideal, ZZ) := (I,n) -> weakLin(module I, n)

weakLin(Module, ZZ,ZZ) := (M, n, d) ->(
        F := res(M, LengthLimit => n);
    <<netList (for i from 1 to n list 
	apply(d, j-> (trim minors(j+1, F.dd_i))))<<endl;
	    )
weakLinList = method()
weakLinList (Module, ZZ) := (M, n) ->(
    S = ring M;
    F := res(M, LengthLimit => n);
    I := trim ideal F.dd_1;
    Jtot := I;    
    for i from 1 to n do(
	I = trim ideal F.dd_(i+1);
	t:= (gens I % Jtot == 0);
	<<(i,t)<<endl;
    	Jtot = trim(Jtot+I))
    )


S3 = ZZ/101[a,b,c]
m3 = ideal gens S3

 Examples = hashTable{(1, ideal"a4,b4,c4,ab3"),
    (2,ideal"a4,b4,c4,ab3,cb3"),
    (3, ideal"a4,b4,c4,ab3,bc3"),
    (4, ideal"a4,b4,c4,abc"),
    (5, ideal"a4,b4,c4,ab3, b2c2")}
S2 = ZZ/101[a,b] 
m2 = ideal gens S2

socle = method()
socle Module := M -> (
    R := ring M;
    ker (M ** transpose vars R)
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

colsInIdeal = method()
colsInIdeal(Matrix, Ideal) := (M, I) -> (
    select(numcols M, j -> gens ideal (M_{j}) % I == 0)
)

socleSummands(ChainComplex, ZZ) := o -> (C,m) -> (
    R := ring C;
    mm := ideal gens R;
    socR := ideal(0_R):mm;
    cyclemaps := for i from 1 to m list syz C.dd_i;
    
    socles := apply(cyclemaps, i -> submatrix(i, colsInIdeal(i, socR)));
--    sc := apply(cycles, c -> select(numcols c, i -> gens ideal c_{i} % socR == 0));
    apply(socles, ss ->(
	    if o.Verbose == false then
       numcols ss
          else
       flatten degrees source ss
       ))
--ssc -> #ssc)
)

socleSummands ChainComplex := o -> C ->  socleSummands(C, length C, o)

socleSummands(Module,ZZ) := o -> (M,ell) -> (
    mm := ideal gens ring M;
    ss := socle M;
    s0 := numcols compress (gens ss % (mm*M));
    F := res(M, LengthLimit => ell);
    {s0}|socleSummands (F,o)
    )

hasSocleSummand = method()
hasSocleSummand Module := M -> (
    0 != socleSummands M)

///
restart
loadPackage "SocleSummands"
S = ZZ/101[x,y]
I =ideal"x4,x2y2,y4"
R=S/I
kSi(R,6)
koSi R
socleSummandGeneration R

k=coker vars R
kSS(R,15)
F = res(k, LengthLimit =>5)
F.dd_4
M = image F.dd_3
--(0*M):(ideal gens R)
socleSummands image F.dd_3
socleSummands coker F.dd_4
C = image (F.dd_3)_{0..5}
MSi(C,7)
///

///
restart
S = ZZ/101[x]
I = ideal "x2"
R=S/I
F = res coker vars R
M = image F.dd_2
(0*M):(ideal gens R)
///


--generate a random monomial of given degree
--The distribution is not uniform, find a way to make it so
randomMonomial = method()
randomMonomial(Ring,ZZ) := (R,d) -> (
	n   := #gens R;
	M   := 1_R;
	apply(d, i -> (M = M*(R_(random(n))))); --multiply d random generators
	M
	)


--Generation of ideals
--generates a homogeneous ideal which gives a graded Artin Ring
randomArtinIdeal = method()

randomArtinIdeal(Ring, ZZ, ZZ) := (R,mexp,extragens) -> (
    --maximum exponent is 1+mexp
	n     := #gens R;
	I     := ideal(0_R);
	for i from 0 to n-1 do (I = I + ideal(R_i^(2+random(mexp))));
	for i from 1 to extragens do (I = I + random(2+random(mexp),R));
--	for i from 1 to extragens do (I = I + randomMonomial(R,2+random(mexp)));
	trim I
	)
randomArtinIdeal(Ring, ZZ) := (R,extragens) -> randomArtinIdeal(R,20,extragens)

--generated at same degree
randomArtinDegreeIdeal = method()
randomArtinDegreeIdeal(Ring, ZZ, ZZ) := (R,d,extragens) -> (
	n     := #gens R;
	I     := ideal(0_R);
	for i from 0 to n-1 do (I = I + ideal(R_i^d));
	for i from 1 to extragens do (I = I + randomMonomial(R,d));
	ideal(mingens I)
	)

--generate a random monomial ideal
randomArtinMonomialIdeal = method()
randomArtinMonomialIdeal(Ring, ZZ) := (R,extragens) -> (
	mexp  := 20;
	n     := #gens R;
	I     := ideal(0_R);
	for i from 0 to n-1 do (I = I + ideal(R_i^(2+random(mexp))));
	for i from 1 to extragens do (I = I + randomMonomial(R,2+random(mexp)));
	ideal(mingens I)
	)

randomArtinMonomialIdeal(Ring, ZZ, ZZ) := (R,mexp,extragens) -> (
	n     := #gens R;
	I     := ideal(0_R);
	for i from 0 to n-1 do (I = I + ideal(R_i^(2+random(mexp))));
	for i from 1 to extragens do (I = I + randomMonomial(R,2+random(mexp)));
	ideal(mingens I)
	)

--Should not be called homogeneous, generated at same degree
randomArtinDegreeMonomialIdeal = method()
randomArtinDegreeMonomialIdeal(Ring, ZZ, ZZ) := (R,d,extragens) -> (
	n     := #gens R;
	I     := ideal(0_R);
	for i from 0 to n-1 do (I = I + ideal(R_i^d));
	for i from 1 to extragens do (I = I + randomMonomial(R,d));
	ideal(mingens I)
	)


ran = method(Options => {Check =>true})
ran(ZZ, Ideal) := Ideal => o -> (d,I) ->  ideal (gens I * random(source gens I, (ring I)^{-d}))
ran(ZZ, Ideal, Ideal) := Ideal => o -> (d,I,J) ->  (
    I' := ideal(gens I % J);
    ran(d, I')
    )

ran(List, Ideal) := Ideal => o-> (L,I) ->(
--    print o.Check;
    L' := sort L;
    g := ran (L'_0,I);
    apply(#L'-1, i -> g = ran(L'_(i+1), I,g) + g);
    if (o.Check == true and codim g<#L) then return "failed" else return g)

ran(ZZ,Ring) := ideal => o-> (d,S) ->(
    ideal random(S^1, S^{-d}))

ran(List, Ring) := ideal => o-> (L,S) ->(
    ideal random(S^1, S^(-L))
)
--Generation of Artin Rings
--just a homogeneous Artin ring
randomArtinRing = method()
randomArtinRing(Ring, ZZ, ZZ) := (R,n,extragens) -> (
	S := R[vars(0..(n-1))];
	I := randomArtinIdeal(S,extragens);
	S/I
	)

--mexp version
randomArtinRing(Ring, ZZ, ZZ, ZZ) := (R,n,mexp,extragens) -> (
	S := R[vars(0..(n-1))];
	I := randomArtinIdeal(S,mexp,extragens);
	S/I
	)

--generated in a single degree
randomArtinDegreeRing = method()
randomArtinDegreeRing(Ring, ZZ, ZZ, ZZ) := (R,n,d,extragens) -> (
	S := R[vars(0..(n-1))];
	I := randomArtinDegreeIdeal(S,d,extragens);
	S/I
	)

--Rings based on Monomial Ideals
randomArtinMonomialRing = method()
randomArtinMonomialRing(Ring, ZZ, ZZ) := (R,n,extragens) -> (
	S := R[vars(0..(n-1))];
	I := randomArtinMonomialIdeal(S,extragens);
	S/I
	)

--mexp version
randomArtinMonomialRing(Ring, ZZ, ZZ, ZZ) := (R,n,mexp,extragens) -> (
	S := R[vars(0..(n-1))];
	I := randomArtinMonomialIdeal(S,mexp,extragens);
	S/I
	)

--Rings based on monomial ideals generated in a single degree
randomArtinDegreeMonomialRing = method()
randomArtinDegreeMonomialRing(Ring, ZZ, ZZ, ZZ) := (R,n,d,extragens) -> (
	S := R[vars(0..(n-1))];
	I := randomArtinDegreeMonomialIdeal(S,d,extragens);
	S/I
	)


///
restart
load "SocleSummands.m2"
S = ZZ/101[a,b,c,d]
I = ideal(a^3,b^3,c^3,d^3,a*b*d,c*d^2,b^2*d,a^2*d,a*c^2,a*d^2,a^2*c,c^2*d,b*c*d)
R = S/I
F = res(coker vars R)
time J = I:ideal vars S
time Rbar = R/sub(J, R)
time L = apply(length F, i->F.dd_i);
time Lbar = L/(phi ->phi**Rbar);
    apply(#Lbar, i ->(
	    s = syz(Lbar_i, DegreeLimit => max flatten degrees source Lbar_i);
	    numcols s ===numcols (L_i)
	    ))
s = syz (Lbar_3);
numcols s
i= 3
i
///


hasSocleSummand1 = method(Options => {Numinfo => false})
hasSocleSummand1(Module) := o-> M -> (
    --this version REQIRES M to be a submodule of a free module ambient M that shares the same socle!
    --this is the case when M is a module of cycles in a minimal free resolution.
	R := ring M;
	mR := ideal vars R;
	socR := (ideal 0_R):mR;
	socSummands := gens(socR*ambient M) % (mR*M);
        if o.Numinfo then numgens prune coker socSummands 
        else socSummands !=0
    	)

///
restart
uninstallPackage"SocleSummands"
loadPackage("SocleSummands", Reload=>true)

S = ZZ/101[a..d]
I = monomialIdeal(a^3,a^2*b,a*b^2,b^3,a^2*c,a*b*c,c^3,a^2*d,d^3)
R = S/I
F = res coker vars R
M = image F.dd_2
hasSocleSummand M
hasSocleSummand1 (M, Numinfo =>true)
kSi(F, Numinfo =>true)
kSi(I, Numinfo =>true)
///
findSocleSummands1 = method(Options =>{Numinfo=>false})
findSocleSummands1(Module) := o-> M -> (
    --this version REQIRES M to be a submodule of a free module ambient M that shares the same socle!
    --this is the case when M is a module of cycles in a minimal free resolution.
	R := ring M;
	mR := ideal vars R;
   	N := (0_R*M) : mR; 
	(gens N) % (mR*M)
)

///
restart
loadPackage("SocleSummands", Reload => true)
S = ZZ/32003[a,b,c]
R = S/ideal(a^4,b^4,c^4,a*b^3,a^2*c^2)
k = coker vars R
use S
R = S/ideal"a2,b2,c2"
kSi(R,7, Numinfo => true)
kSS(R,7)
F = res(k, LengthLimit => 12)
time hasSocleSummand1(image F.dd_4)
time hasSocleSummand(image F.dd_4)
M = image F.dd_9;
time hasSocleSummand M
kSi R
kSi(R,10)

--
kk = ZZ/101
S = kk[a,b,c]
ell = ideal(a^3,a^2*b,a*b^2,b^3,a^2*c,a*b*c)
R = S/ell
K = koszul vars R
M = ker K.dd_1
M' = prune M
hasSocleSummand M'
hasSocleSummand M
hasSocleSummand1 M
///

kSi = method(Options=>{Numinfo => false})
kSi ChainComplex := o -> F -> (
    --does not assume that F is a resolution, eg Koszul complex
    ell := length F;
    apply(ell-1, i-> (
	    hasSocleSummand1(image syz F.dd_(i+1),o)
	       ))
       )

kSRes = method(Options=>{Numinfo => false})
--assumes that F.dd_(i+1) = syz F.dd_i.
--does not test coker F.dd_1; starts with image F.dd_2
kSRes ChainComplex := o -> F -> (
    ell := length F;
    apply(ell-1, i-> (
	    hasSocleSummand1(image F.dd_(i+2),o)
	       ))
       )
   
kSRes(ChainComplex, ZZ) := o-> (F,maxlength) -> (
    --this version controls the number of terms to compute
    --even when the resolution F is too long or too short.
    if length F < maxlength then F' := res(coker F.dd_1, 
	                             LengthLimit =>maxlength)
    else F' = F;
    apply(maxlength-2, i-> (
	    hasSocleSummand1(image F'.dd_(i+2),o)
	       ))
       )

kSi(Ideal, ZZ) := o -> (I, maxres) -> (
    R := ring I/I;
    kSRes(res(coker vars R, LengthLimit => maxres),o)
    )

kSi Ideal := o -> I -> kSi (ring I/I, o)
    
kSi(Module, ZZ) := o-> (M, maxres) ->(
	F := res(M, 
		LengthLimit => maxres + 1);
	kSRes(F,o))

kSi(QuotientRing, ZZ) := o -> (R, maxres) -> 
    kSi(coker vars R, maxres, o)

kSi QuotientRing := o-> R -> (
    maxres := #gens R + 2;
    F := res(coker vars R, LengthLimit => maxres);
    kSRes(F,o))

///
restart
loadPackage "SocleSummands"
S = ZZ/101[a..d]
I = monomialIdeal(a^3,a^2*b,a*b^2,b^3,a^2*c,a*b*c,c^3,a^2*d,d^3)
R = S/I
kSi(I,6,Numinfo => true)
kSi(S/I,5)
kSi (S/I)
F= res(coker vars R, LengthLimit =>7)
elapsedTime kSi F -- 9.3 sec
elapsedTime kSRes F --7.1 sec
///


///
restart
load "SocleSummands.m2"
needsPackage "DGAlgebras"
S = ZZ/101[a,b,c]
 I1 = ideal (a^4,b^4,c^4, a*b^3)
koSi (S/I1)
kSi (S/I1)
isGolod (S/I1)
--this example shows that in a non-golod ring the 
--socle summands in the Koszul complex can become
--submodules in the max ideal times the syzygy in
--the minimal free res.

 I2= ideal (a^4,b^4,c^4, a*b^3,c*b^3) 
 kSi (I2, 4)
use S
 kkSi(I2)
 I3=ideal (a^4,b^4,c^4, a*b^3, b*c^3)
 I4= ideal (a^4,b^4,c^4,a*b*c)
 I5 = ideal (a^4,b^4,c^4, a*b^3, b^2*c^2)

kSi(QuotientRing) := o-> R -> (kSi(R, 4 + #gens R)
    )
///

--determines if a is in the numerical semigroup gen by L
--elements of L are assumemed to be postive integers
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

--The numerical semigroup of kS_i indexes
kSS = method()
kSS(QuotientRing, ZZ) := (R,n) -> (
    	mm := ideal gens ambient R;
	I := ideal R;
	if 0 != gens(mm*(I:mm)) % (mm*I) then return {2,3};	
	K := kSi(R,n); --compute kSi
	L := {}; --empty list
	newgen := null;
	for i from 1 to (n-1) do (
		newgen = K_i and (not inSemigroup(i,L));
		if newgen then (L = L | {i});
		);
	if (L == {}) then (L = {0});
	L
	)
kSS(Ideal, ZZ) := (I,n) -> kSS((ring I)/I, n)

--The k index. this is only an approximation
--if n is too small then the output is infinity
kSI = method()
kSI(QuotientRing, ZZ) := (R,n) -> (
	K := kSi(R,n);
	L := infinity;
	for i from 1 to (n-1) do (
		if K_i then (
			L = i;
			break;
			);
		);
	L
	)

///
viewHelp MSi
///
--Testing the (MSi) Property
--Less useful I think
MSi = method(Options => {Numinfo => false})
MSi(Module, ZZ) := o ->(M, maxres) -> (
	F := res(M, LengthLimit => maxres + 1);
	<<betti F<<endl;
	scan(maxres, i -> (
		if i == 0 then 
		   <<hasSocleSummand(M)<<endl --maybe slow but is now correct.
		else (
		    N := image F.dd_i;
		    b := hasSocleSummand1 N;
		    if b == false then <<false<<endl
		        else <<hasSocleSummand N<<endl)))
	)

MSiList = method()
MSiList(Module, ZZ) := (M, maxres) -> (
	F := res(M, LengthLimit => maxres + 1);
	apply(maxres, i -> (
		if i == 0 then false
		else (
		    N := image F.dd_i;
		    b := hasSocleSummand1 N;
		    if b == false then false
		    else hasSocleSummand N 
		     )
	)))
///
MSiList(R^1/Ilist_4, 6)
///
--Compute the induced map Tor_i(mI,k) -> Tor_i(I,k), R = S/I
TorMap = method()
TorMap(QuotientRing,ZZ) := (R,i) -> (
	I := ideal R;
	S := ambient R;
	m := ideal vars S;
	f := inducedMap(image gens I, image mingens (m*I));
	k := coker vars S;
	F := res(f, LengthLimit => (numgens S) + 2);
	G := F**k;
	M := (ker (source G).dd_(i))/(image (source G).dd_(i+1));
	N := (ker (target G).dd_(i))/(image (target G).dd_(i+1));
	b := inducedMap(N,M,G_i)
	)
--the following prints a list of booleans: true in the i-th spot means
--Tor_(i+1)(mI,k) ->Tor_(i+1)(I,k) is surjective, where i = 0,1... 
--The list has n-1 entries.
--Conjecture: if true at the j-th place in the output (counting from 0!), 
--then R is not kS_(n-j).
TorMapk = method()
TorMapk QuotientRing := List => R -> (
	I := ideal R;
	S := ambient R;
	m := ideal vars S;
	f := inducedMap(image gens I, image mingens (m*I));
	k := S/ideal(vars S);
	F := res(f, LengthLimit => (numgens S) + 2);
    	apply(length target F, i-> rank(k**F_(i+1)) == rank target (F_(i+1))
	    )
	    )
///
kSS(R,7)
///
TorMapConjecture = method()
TorMapConjecture QuotientRing := Boolean => R ->(
    n := numgens R;
    L := TorMapk R;
    PL := positions(L, ell-> ell == true);
    K := kSi(R,n+1);
-*
    PK := positions(K, k-> k == true);
    all(PK, i-> all(PL,ell-> ell !=n-i))
*-
    P := positions(L, b -> b==true);
    Pdual :=  apply(P,i->n-i);
    all(Pdual, i-> K_i == false)
    )
    
///
--counterexample to the Tor map conjecture, 8/17/2020:
S = kk[a,b,c,d]
i = ideal(a^3,b^3,c^3,d^3,a*b*c,b^2*d,a^2*b)
R = S/i
TorMapk(S/i)
///
    
///
restart
loadPackage("SocleSummands", Reload =>true)
S = ZZ/101[a,b,c]
I = {ideal"a4,b4,c4,abc", ideal(a^2,b^2,c^2),ideal"a4,b4,c4,ab3,bc3",
    (ideal"a2,b2,c2"*ideal vars S),ideal"a4,b4,c4,a2b2,a2c2"}
apply(I,i->TorMapConjecture(S/I))


R = S/I
TorMapConjecture R
use S
TorMapk(S/ideal(a^2,b^2,c^2))
TorMapk(R)
TorMapConjecture R


I = ideal"a4,b4,c4,ab3,bc3"
R = S/I
TorMapk(R)

kSS(R,7)
TorMap(R,1)
use S
m = ideal vars S
f = inducedMap(image gens I, image mingens (m*I))
F = res f

use S
R = S/ideal"a3,b3,c3,a2b2c2"
TorMapk R
kSS(R, 7)
TorMapConjecture R
kSS(R,7)
use S
R = S/(ideal"a2,b2,c2"*ideal vars S)
TorMapk R
TorMapConjecture R
kSS R

restart
loadPackage("SocleSummands", Reload =>true)
S = ZZ/101[a,b,c]
m = ideal vars S
I = ideal"a4,b4,c4,a2b2,a2c2"
R = S/I
betti res (S^1/I)
betti res (S^1/(m*I))

TorMap(R,2)
TorMapk R

///
--Compute the induced map Tor_i(M,P) -> Tor_i(N,P) given M -> N
TorMap(Matrix, Module, ZZ) := (f,P,i) -> (
	F := res(f, LengthLimit => (i+1));
	G := F**P;
	M := (ker (source G).dd_(i))/(image (source G).dd_(i+1));
	N := (ker (target G).dd_(i))/(image (target G).dd_(i+1));
	b := inducedMap(N,M,G_i)
	)
--check for which i the map of resolutions of mmI to I is surjective.
--TorMap(QuotientRing, ZZ) 

--Compute the maps above and return a list of such maps

TorMaps = method()
TorMaps(Matrix, Module, ZZ) := (f, P, iend) -> (
	M := source f;
	N := target f;
	C := res(M, LengthLimit => (iend + 1));
	D := res(N, LengthLimit => (iend + 1));
	F := extend(D,C, matrix f);
	G := F ** P;
	L := {};
	for i from 0 to iend do (
		A := (ker (source G).dd_(i))/(image (source G).dd_(i+1));
		B := (ker (target G).dd_(i))/(image (target G).dd_(i+1));
		c := inducedMap(B,A,G_i);
		L = L|{c};
		);
	L
	)

--Checks the surjectivity of the Tor maps above
isTorSurjective = method()
isTorSurjective(QuotientRing) := R -> (
	I := ideal R;
	S := ambient R;
	m := ideal vars S;
	f := inducedMap(image gens I, image mingens (m*I));
	k := coker vars S;
	F := TorMaps(f,k,(#gens R)+1);
	F/isSurjective
	)

isTorSurjective(QuotientRing, ZZ) := (R,iend) -> (
	I := ideal R;
	S := ambient R;
	m := ideal vars S;
	f := inducedMap(image gens I, image mingens (m*I));
	k := coker vars S;
	F := TorMaps(f,k,iend+1);
	F/isSurjective
	)


--given a strictly increasing list this method returns the list that increases
--the entry of the smallest index by one for which this increase still results
--in a strictly increasing list
-*
def consecutive(L):
	c = True
	for i in range(len(L) - 1):
		if L[i+1]-L[i] != 1:
			c = False
			break
	return c
*-
isConsecutive = method()
isConsecutive(List) := L -> (
	val := true;
	for i from 0 to (#L-2) do (
		val = val and ( L#(i+1)-L#i == 1 );
		);
	val
	)


nextList = method()
nextList(List) := (L) -> (
	M := {};
	if isConsecutive(L) then (
		M = toList(0..(#L-2)) | {L#(-1) + 1};
		)
	else (
		M = nextList(L_{0..(#L-2)});
		M = M | {L#(-1)};
		);
	M
	)


--given a list of one less entry than the number of generators of the
--ring R this method returns the monomial corresponding to the
--stars and bars representation of the monomial
ListToMultiDegree = method()
ListToMultiDegree(Ring, List, ZZ) := (R, L, d) -> (
	deg := {L#0};
	for i from 1 to (#L -1) do (
		deg = deg | {L#(i) - L#(i-1) - 1};
		);
	deg = deg | {d - L#(-1)};
	R_deg
	)


--returns a list of the degree d monomials of R that are
--not of the form x_i^d
monomialList = method()
monomialList(Ring, ZZ) := (R, d) -> (
	n        := #(gens R);
	monoList := {};
	N        := d+n-1; --stars and bars total
	degList  := toList(0..(n-2)); --first multidegree
	mono     := null; --monomial
	while (degList#(-1) < N) do (
		mono     = ListToMultiDegree(R, degList, N-1);
		monoList = monoList | {mono};
		degList  = nextList(degList);
		);
	for i from 0 to (n-1) do (
		monoList = delete(R_i^d, monoList);
		);
	monoList
	)


--returns a list of all monomial ideals of the form
--(x_1^d,...,x_n^d, extragens more monomials of deg d)
monomialIdealList = method()
monomialIdealList(Ring, ZZ, ZZ) := (R, d, extragens) -> (
	n         := #(gens R);
	monos     := monomialList(R,d);
	standard  := {};
	idealList := {};
	for i from 0 to (n-1) do (
		standard = standard | {R_i^d};
		);
	combo := toList(0..(extragens - 1));
	while combo#(-1) < #monos do (
		idealList = idealList | {ideal(standard | monos_combo)};
		combo = nextList(combo);
		);
	idealList
	)

///
S = ZZ/101[a,b,c,d]
L = monomialIdealList(S,3,3)
#L
apply(L, i-> (
	t = TorMapConjecture(S/i);
	if t === false then <<i<<endl<<flush;
	)
)
time all(apply(L, i-> TorMapConjecture(S/i)), t -> t===true)

i = ideal"a3,b3,c3,d3,abc,b2d,a2b"
TorMapConjecture (S/i)
kSS(S/i, 7)
TorMapk(S/i)

///


--determines if list of extra generators is minimial under symmetric
--group action
isReducedList = method()
isReducedList(List) := L -> (
	L     = sort L; --sort the list
	M    := L; --dummy value
	S    := ring(L#0); --ambient ring
	n    := #(gens S);
	Perm := permutations n; --permutations
	flag := true; --minimality status, innocent until proven guilty
	f    := null; --placeholder for map in loop
	for g in Perm do (
		f = map(S,S,apply(n,i -> S_(g#i))); --permuation ring map
		M = apply(L,i->f(i)); --apply f to L
		M = sort M; --sort for comparison
		if M < L then (
			flag = false;
			break;
			); --if L is not minimal then break
		);
	flag --return the minimality status
	)
S = ZZ/32003[x,y,z]
lin7 = {monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^2*y^4*z,y^6*z,x
      ----------------------------------------------------------------------------------------------------
      ^5*z^2,x^4*y*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z
      ----------------------------------------------------------------------------------------------------
      ^5,x*z^6,y*z^6,z^7), monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z
      ----------------------------------------------------------------------------------------------------
      ,x^4*y^2*z,x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y
      ----------------------------------------------------------------------------------------------------
      ^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^
      ----------------------------------------------------------------------------------------------------
      4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*
      ----------------------------------------------------------------------------------------------------
      y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^
      ----------------------------------------------------------------------------------------------------
      4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,y^
      ----------------------------------------------------------------------------------------------------
      3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,
      ----------------------------------------------------------------------------------------------------
      y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^3*y^3*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,
      ----------------------------------------------------------------------------------------------------
      x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^3*y^3*z,x^2*y^4*
      ----------------------------------------------------------------------------------------------------
      z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,
      ----------------------------------------------------------------------------------------------------
      x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^2*y^4*z,y^6*z,x^
      ----------------------------------------------------------------------------------------------------
      5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*z^3,x^3*z^4,x^2*y*
      ----------------------------------------------------------------------------------------------------
      z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*
      ----------------------------------------------------------------------------------------------------
      y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,
      ----------------------------------------------------------------------------------------------------
      x^3*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*
      ----------------------------------------------------------------------------------------------------
      y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,
      ----------------------------------------------------------------------------------------------------
      x^3*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x*
      ----------------------------------------------------------------------------------------------------
      y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*z^3,x^3*z^
      ----------------------------------------------------------------------------------------------------
      4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^3,x^3*z^
      ----------------------------------------------------------------------------------------------------
      4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^3,x^3*z^
      ----------------------------------------------------------------------------------------------------
      4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,
      ----------------------------------------------------------------------------------------------------
      x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^3*y^3*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,
      ----------------------------------------------------------------------------------------------------
      x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^3*y^3*z,x^2*y^4*
      ----------------------------------------------------------------------------------------------------
      z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*z^3,x^3*
      ----------------------------------------------------------------------------------------------------
      z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^3*y^3*z,x^2*y^4*
      ----------------------------------------------------------------------------------------------------
      z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^
      ----------------------------------------------------------------------------------------------------
      4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,y^4*z^3,x^
      ----------------------------------------------------------------------------------------------------
      3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*
      ----------------------------------------------------------------------------------------------------
      y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,
      ----------------------------------------------------------------------------------------------------
      x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^
      ----------------------------------------------------------------------------------------------------
      2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^3,
      ----------------------------------------------------------------------------------------------------
      x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^3*y^3*z,x^2*y^4*
      ----------------------------------------------------------------------------------------------------
      z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^
      ----------------------------------------------------------------------------------------------------
      4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^3*y^3*z,x^2*y^4*
      ----------------------------------------------------------------------------------------------------
      z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^
      ----------------------------------------------------------------------------------------------------
      3,x*y^3*z^3,y^4*z^3,x^3*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,y^4*z^3,x^
      ----------------------------------------------------------------------------------------------------
      3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^3,y^4*z^3,
      ----------------------------------------------------------------------------------------------------
      x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^
      ----------------------------------------------------------------------------------------------------
      2*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^
      ----------------------------------------------------------------------------------------------------
      2*z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^
      ----------------------------------------------------------------------------------------------------
      3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^
      ----------------------------------------------------------------------------------------------------
      3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^
      ----------------------------------------------------------------------------------------------------
      3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^
      ----------------------------------------------------------------------------------------------------
      3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^
      ----------------------------------------------------------------------------------------------------
      3*z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,y^4*
      ----------------------------------------------------------------------------------------------------
      z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^2*y^4*z,
      ----------------------------------------------------------------------------------------------------
      x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,y^4*z^
      ----------------------------------------------------------------------------------------------------
      3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^4*y^2*z,x^3*y^3*z,x^2*y^4*
      ----------------------------------------------------------------------------------------------------
      z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*
      ----------------------------------------------------------------------------------------------------
      y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^
      ----------------------------------------------------------------------------------------------------
      3,x^2*y^2*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^
      ----------------------------------------------------------------------------------------------------
      3,x^2*y^2*z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^
      ----------------------------------------------------------------------------------------------------
      3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*
      ----------------------------------------------------------------------------------------------------
      z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*
      ----------------------------------------------------------------------------------------------------
      z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*
      ----------------------------------------------------------------------------------------------------
      z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^
      ----------------------------------------------------------------------------------------------------
      3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^
      ----------------------------------------------------------------------------------------------------
      3,x*y^3*z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^3*z^
      ----------------------------------------------------------------------------------------------------
      3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^3*y^2*z^2,x^2*y^3*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^2*z^3,x*y^3*
      ----------------------------------------------------------------------------------------------------
      z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^
      ----------------------------------------------------------------------------------------------------
      2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^3,x^2*y^
      ----------------------------------------------------------------------------------------------------
      2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^2*y^2*z^3,x*y^
      ----------------------------------------------------------------------------------------------------
      3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^
      ----------------------------------------------------------------------------------------------------
      3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,y^2*z^5,x*z^6,y*z^6,z^7
      ----------------------------------------------------------------------------------------------------
      ), monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3
      ----------------------------------------------------------------------------------------------------
      *z,x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y
      ----------------------------------------------------------------------------------------------------
      *z^3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7
      ----------------------------------------------------------------------------------------------------
      ), monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3
      ----------------------------------------------------------------------------------------------------
      *z,x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y
      ----------------------------------------------------------------------------------------------------
      *z^3,x^2*y^2*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7
      ----------------------------------------------------------------------------------------------------
      ), monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3
      ----------------------------------------------------------------------------------------------------
      *z,x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y
      ----------------------------------------------------------------------------------------------------
      *z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*z^6,z^7),
      ----------------------------------------------------------------------------------------------------
      monomialIdeal(x^7,x^6*y,x^5*y^2,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^6*z,x^5*y*z,x^4*y^2*z,x^3*y^3*z,
      ----------------------------------------------------------------------------------------------------
      x^2*y^4*z,x*y^5*z,y^6*z,x^5*z^2,x^4*y*z^2,x^3*y^2*z^2,x^2*y^3*z^2,x*y^4*z^2,y^5*z^2,x^4*z^3,x^3*y*z^
      ----------------------------------------------------------------------------------------------------
      3,x^2*y^2*z^3,x*y^3*z^3,y^4*z^3,x^3*z^4,x^2*y*z^4,x*y^2*z^4,y^3*z^4,x^2*z^5,x*y*z^5,y^2*z^5,x*z^6,y*
      ----------------------------------------------------------------------------------------------------
      z^6,z^7)};

///
restart
debug loadPackage "SocleSummands"
numgens lin7_0
numcols basis (7,S)
///

///
restart
loadPackage "SocleSummands"
needsPackage "Points"
S = ZZ/101[a,b]
mm = ideal gens S
R = S/(mm^[2])^2
mmR = ideal gens R
weakLin(module(mmR^2), 4)
weakLin(mmR^2,4)
R = ringOfRandomPoints(2,11)
mmR = ideal gens R
weakLin (mmR^3, 4)
MSi(module(mmR^2), 7)

test = method(Options => {Bound => 7})
test ZZ := o-> n  ->(
    print o.Bound;
)
test(1, Bound =>5)

///

ringOfRandomPoints  = (r,n)-> (
    --artinian reduction of  the ring of n random points in PP^r
    needsPackage "Points";
    x:= symbol x;
    S := ZZ/32003[x_0..x_r];
    I := points randomPointsMat(S,n);
    R := coefficientRing S[(gens S)_{1..r}];
    phi := map(R,S,{random(1,R), R_0..R_(r-1)});
    assert(codim phi I == r);
    R/phi I)


resOfRandomPoints = (r,n) -> (
    Rbar := ringOfRandomPoints(r,n);
    S = ambient Rbar;
    use S;
    << minimalBetti (S^1/ideal Rbar)<<", "<<betti (F := res(coker vars Rbar, LengthLimit => r+1))<<endl;
    F
    )
linearFormsInRes = I ->(
    mm := ideal vars ring I;
    F := res module I;
    apply(length F , i-> t :=(0 != compress ((gens ideal F.dd_(i+1))%mm^2)))
)
      

beginDocumentation()
--Documentation
-*
restart
loadPackage("SocleSummands", Reload=> true)
installPackage "SocleSummands"
check SocleSummands
Todo list:
	"randomMonomial",
	"randomArtinIdeal",
	"randomArtinDegreeIdeal",
	"randomArtinMonomialIdeal",
	"randomArtinDegreeMonomialIdeal",
	"randomArtinRing",
	"randomArtinDegreeRing",
	"randomArtinMonomialRing",
	"randomArtinDegreeMonomialRing",
*-
doc ///
	Key
		hasSocleSummand
		(hasSocleSummand, Module)
	Headline
		Determines if residue field is a direct summand of input
	Usage
		 hasSocleSummand M
	Inputs
		M:Module
			A module over an Artin local ring
	Outputs
		:Boolean
			Whether the residue field is a direct summand of M
	Description
		Text
			Takes in a module over a local Artin ring and determines
			is the residule field of the ring is a direct summand of
			M
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			M = coker matrix {{a*b,b*c},{c,a*b*c}}
			hasSocleSummand M

///

doc///
Key
 socle
 (socle, Module)
Headline
 compute the socle
Usage
 s = socle M
Inputs
 M:Module
Outputs
 s:Module
   submodule of M
Description
  Text
   returns  ker (M ** transpose vars R).
  Example
   R = ZZ/101[a,b]/(ideal"a3,a2b,b4")
   socle R^1
SeeAlso
 socleSummands
///

doc ///
Key
 socleSummands
 (socleSummands, Module)
 (socleSummands, Module, ZZ)
 (socleSummands, ChainComplex, ZZ) 
 (socleSummands, ChainComplex)
Headline
 counts the socle summands of a module or cycles in a complex
Usage
 s =  socleSummands M
 s =  socleSummands(M,ell)
 ss =  socleSummands(C,ell)
 ss =  socleSummands C
Inputs
 M : Module
 C: ChainComplex
 ell: ZZ
Outputs
 s: ZZ
 ss: List
  of ZZ
Description
  Text
   socleSummands M returns the number of generators of the socle of a module M.
   
   socleSummands(M,ell) returns the number of socle summands a module M and the first ell cycle modules in 
   a resolution of M.
   
   socleSummand(C,ell)  returns the number of socle summands in the first ell cycle modules of a complex C.
   
   socleSummands C  = socleSummands(C, length C).
  Example
   R = ZZ/101[a,b,c]/(ideal(a,b,c))^2
   socleSummands(coker vars R, 5)
   socleSummands koszul vars R
   socleSummands (res (coker vars R, LengthLimit => 5))
SeeAlso
 socle
///

doc ///
	Key
		randomArtinIdeal
		(randomArtinIdeal, Ring, ZZ, ZZ)
		(randomArtinIdeal, Ring, ZZ)
	Headline
		Produces a random Artinian ideal
	Usage
		randomArtinIdeal(R, m, n)
		randomArtinIdeal(R, n)
	Inputs
		R:Ring
		m:ZZ
			that bounds the degree of the generators
		n:ZZ
			that controls the number of minimal generators
	Outputs
		:Ideal
		    with a minimal set of generators
		    
	Description
		Text
			Let r be the number of generators of R. The function 
			constructs a random homogenous Artinian ideal in R with 
			at most r + n generators. The degree of each generator 
			is at most m + 1. If no value of m is given, the 
			default value is 20.
			
			Note that the generators of the ideal includes random powers of 
			the generators of R.
		Example
			S = ZZ/32003[x,y,z,w]
			I = randomArtinIdeal(S,6,4)
			numgens I
			L = flatten entries mingens I;
			apply(L, v -> first degree v)
    	SeeAlso
	     randomMonomial
	     randomArtinDegreeIdeal
	     randomArtinMonomialIdeal
	     randomArtinDegreeMonomialIdeal
	     randomArtinRing
	     randomArtinDegreeRing
	     randomArtinMonomialRing
	     randomArtinDegreeMonomialRing
    
///

doc ///
	Key
		kSi
		(kSi, QuotientRing)
		(kSi, QuotientRing, ZZ)
	Headline
		A boolean list of when k is a direct summand of its syzygy
	Usage
		kSi R
		kSi(R, n)
	Inputs
		R:QuotientRing
			An Atrin local ring
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:List
			List of booleans when k is a direct summand of its
			syzygy
	Description
		Text
			Takes in an Artin local ring, R, and computes a free
			resolution of its residue field, k, and determines when
			k is a direct summand of its ith syzygy. Outputs a list
			of boolean values for when it is.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			kSi(R,6)
///


doc ///
	Key
		kSS
		(kSS, QuotientRing, ZZ)
	Headline
		A minimial list of generators for kSS(R)
	Usage
		kSS(R, n)
	Inputs
		R:QuotientRing
			An Atrin local ring
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:List
			List of generators for the numerical semigroup generated
			by the indexes for which k is a direct summand of its
			syzygies of the corresponding index
	Description
		Text
			Takes in an Artin local ring, R, and computes a free
			resolution of its residue field, k, and determines when
			k is a direct summand of its ith syzygy. For the idexes
			that return true we return the sublist of minimal
			generators for the numerical semigroup they generate
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			kSS(R, 7)
///


doc ///
	Key
		kSI
		(kSI, QuotientRing, ZZ)
	Headline
		The (kS) index of R
	Usage
		kSI(R, n)
	Inputs
		R:QuotientRing
			An Atrin local ring
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:ZZ
			the infimum of indexes i for which R is (kSi) for i > 0
	Description
		Text
			Takes in an Artin local ring, R, and integer n and
			computes a free resolution of its residue field, k, up
			to index n and determines when k is a direct summand of
			its ith syzygy. If any positive index results in (kSi)
			being true then this will return the minimum of these
			indexes. Otherwise this method returns infinity.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			kSI(R,7)
///


doc ///
	Key
		MSi
		(MSi, Module, ZZ)
	Headline
		A boolean list of when k is a direct summand of the syzygies of M
	Usage
		MSi(M, n)
	Inputs
		M:Module
			A module over an Atrin local ring, R
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:List
			List of booleans when k is a direct summand of the
			syzygies of M
	Description
		Text
			Takes in a module over an Artin local ring, R, and
			computes a free resolution of its residue field, k, and
			determines when k is a direct summand of its ith
			syzygy. Outputs a list of boolean values for when it is.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			M = coker matrix {{a*b,b*c},{c,a*b*c}}
			MSi(M, 6)
///

doc ///
	Key
		TorMap
		(TorMap, QuotientRing, ZZ)
		(TorMap, Matrix, Module, ZZ)
	Headline
		Computes the induced map on Tor_i
	Usage
		TorMap(R, i)
		TorMap(f, P, i)
	Inputs
		R:QuotientRing
			A module over an Atrin local ring, R
		f:Matrix
			Module map from M to N
		P:Module
			Fuctor applied is Tor_i(-,P)
		i:ZZ
			Index of Tor
	Outputs
		:Matrix
			The induced map Tor_i(M,P) -> Tor_i(N,P)
	Description
		Text
			Given an Artin local ring R with ambient ring S
			and defining ideal I in S and integer i, TorMap(R,i)
			computes the induced map Tor_i(mI,k) -> Tor_i(I,k)
			where m is the maximal ideal generated by the variables
			of S and k is the residue field of R thought of as an
			S-module
			
			Given modules M,N,P over an Artin local ring R, an
			integer i, and a matrix f:M->N thought of as a module
			map, TorMap(f,P,i) computes the induced map
			Tor_i(M,P) -> Tor_i(N,P)
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			M = R^3
			N = coker matrix {{a*b},{b*c}}
			P = coker matrix {{a},{b},{c}}
			F = matrix {{a^2,b^2,c^2},{a*b,b*c,c*a}}
			f = inducedMap(N,M,F)
			TorMap(f,P,3)
			TorMap(R,3)
///

doc ///
	Key
		TorMaps
		(TorMaps, Matrix, Module, ZZ)
	Headline
		Computes a list of the induced maps on Tor_i
	Usage
		TorMaps(f, P, iend)
	Inputs
		f:Matrix
			Module map from M to N
		P:Module
			Fuctor applied is Tor_i(-,P)
		iend:ZZ
			Maximal index of Tor
	Outputs
		:List
			List of the induced maps Tor_i(M,P) -> Tor_i(N,P) for
			i from 0 to iend
	Description
		Text
			Given modules M,N,P over a ring R, an integer i, and a
			module map f:M->N, TorMaps(f,P,i) computes a list of
			the induced maps Tor_i(M,P) -> Tor_i(N,P) for i from 0 
			to iend
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			M = R^3
			N = coker matrix {{a*b},{b*c}}
			P = coker matrix {{a},{b},{c}}
			F = matrix {{a^2,b^2,c^2},{a*b,b*c,c*a}}
			f = inducedMap(N,M,F)
			TorMaps(f,P,3)
///


doc ///
	Key
		isTorSurjective
		(isTorSurjective, QuotientRing)
		(isTorSurjective, QuotientRing, ZZ)
	Headline
		List surjectivity status of the induce Tor maps
	Usage
		isTorSurjective(R)
		isTorSurjective(R, n)
	Inputs
		R:QuotientRing
			Specified ring
		n:ZZ
			Maximal index for which surjectivity status will be
			computed
	Outputs
		:List
			List of the surjectivity status of the induced maps 
			Tor_i(M,P) -> Tor_i(N,P) for i from 0 to n
	Description
		Text
			For a local Artin ring R with ambient polynomial ring S
			and defining ideal I in S. Let m be the maximal ideal of
			S generated by the variables, the preimage of the
			maximal ideal of R. As S-modules we have the inclusion
			map mI -> I which induces maps Tor_i(M,P) -> Tor_i(N,P).
			This method computes a list booleans for i from 0 to n
			determining if the induced map is surjective.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			isTorSurjective(R,3)
///

test1 = method()
test1(Ring, ZZ, ZZ, ZZ) := (k,n,d,extragens) -> (
	S         := k[vars(0..(n-1))]; --Polynomial Ring
	monos     := monomialList(S,d);
	standard  := {};
	dataList  := {};
	I         := null; --blank ideal
	J         := null; --blank ideal	
	K         := null; --blank list
	R         := null;
	Rbar      := null;
	counter   := 0; --counter for debugging purposes
	for i from 0 to (n-1) do (
		standard = standard | {S_i^d};
		);
	combo := toList(0..(extragens - 1));
	while combo#(-1) < #monos do (
		counter = counter + 1; --increment
		--print counter; --display count
		if isReducedList(monos_combo)and counter==1 then (
			print counter;
			I = ideal(standard | monos_combo); --current ideal
			J = I:ideal vars S;
			R = S/I;
			Rbar = R/sub(J,R);
			F := res(coker vars R, LengthLimit => 10);
			--print("okay1");
			
			K = hasSocleSummand1(Rbar, F); --compute kSS
			--print("okay2");
			dataList = dataList | {K}; --change the second arg in kSS to get loger res
			print I;
			print K;
			);
		combo = nextList(combo);
		);
	tally dataList
	)

golodTests = method(Options =>{Verbose => false})
golodTests Ideal :=  o ->I ->  (
    --by paper of Dao and diStefani, this pair of tests is necessary for Golodness,
    --and sufficient for monomial ideals in 3 vars; not known whether its
    --its sufficient for general homogeneous ideals in 3 vars.
    R := ring I;
    v := gens R;
    V := (permutations v)_{0,3,4}; -- the cyclic permutations.
    condition1 := for w in V list(
	gens((I:w_0)*(I:(ideal(w_1,w_2))))%I == 0
	);

    condition2 := for w in V list(
	gens((I:w_0)*(I:w_1)) % (w_2*(I:(ideal(w_0,w_1)))+I) == 0
	);
    a1 := all(condition1, t->t);
    if a1 == false then
        if o.Verbose == false then return false else
	     <<netList I_*<<endl<<"failed condition 1"<<endl;
    a2 := all(condition2, t->t);
    if a2 == false then  
        if o.Verbose == false then return false else 
	      <<netList I_*<<endl<<"failed condition 2"<<endl;    
    true
)    

///
restart
loadPackage "SocleSummands"
R = ZZ/101[a,b,c]
I = ideal"a2,b3,c4,ab,bc"
golodTests I

restart
uninstallPackage "SocleSummands"
time installPackage "SocleSummands"
--check "SocleSummands"
needsPackage "SocleSummands"
viewHelp SocleSummands
time test1(ZZ/101, 3,3,1)
///

///
dist1 = (e1,e2) -> sum apply(#e1, i-> max (e1_i,e2_i)) - sum e1
--e1,e2 should be lists of ZZ of the same length and same sum d, representing two monomials of degree d.
--the distance is the degree of the lcm minus d.
-- this is not exported.

dist = method()
dist(List, List) := ZZ => (E1,E2) -> (
    --E1,E2 should each be a list of lists of ZZ, representing monomials of the same degree.
    --the distance is the minimal distance between monomials in the set.
    min flatten apply(E1, e1-> apply(E2, e2-> dist1(e1,e2)))
	)
///


kkSi = method(Options=>{Nonminimal => false})
kkSi Ring := o-> R -> (
    n := numgens R;
    K := koszul vars R;
--    F := res(coker vars R, LengthLimit =>n, 
--	                 FastNonminimal => o.Nonminimal);
--    <<kSi(F)<<endl;
    for i from 0 to n+1 list (
	if i == 0 then true else
	if i == 1 then hasSocleSummand R^1 else
	hasSocleSummand1 ker K.dd_(i-1)
	)
    )
kkSi Ideal := o -> I -> kkSi(ring I/I)

KtoR = method(Options => {Resolution => 0})
KtoR (Ring, ZZ) := o-> (R,n) -> (
    --print o#Resolution;
    KK := new MutableList;
    KK#0 = koszulComplexDGA R;
    for i from 1 to n do 
	KK#i = killCycles(KK#(i-1), EndDegree => i);
    K := toList KK;
    if o.Resolution == 0 then K else apply(#K, i -> toComplex(K_i, o.Resolution))
    )

///
restart
loadPackage( "SocleSummands", Reload => true)
S = ZZ/101[a,b,c,d]
R = S/(ideal apply(gens S, x -> x^2))^2
KK = KtoR(R,3, Resolution => 0);
F = last KK
(KK_2).natural
degreeLength ((flattenRing oo)_0)
(gens flattenRing oo)/degree

degreeLength ((KK_0).natural)
///
end--

restart
loadPackage("SocleSummands", Reload=> true)
uninstallPackage "SocleSummands"
restart
installPackage "SocleSummands"
--check "SocleSummands"




--test on monomial ideals gen at same degree
--randomArtinDegreeMonomialRing(R,n,d,extragens)
--Things to test
-*
1) (kSi) on the rings
2) The Tor^S_i surjectivity for n+1-i > 0, i < n+1

--Tally the (kSi) patterns
--verify tor surjectivity conjecture
*-
test = method()
test(Ring, ZZ, ZZ, ZZ) := (k,n,d,extragens) -> (
	Ks := {}; --Store (kSi)
	K  := null;
	Ts := {}; --Store surjectiveity of ith Tor map
	R  := ZZ;
	for i from 1 to 100 do (
		R = randomArtinDegreeMonomialRing(k,n,d,extragens); --random ring
		K  = kSi(R);
		Ks = Ks | {K};  --store kSi
		Ts = Ts | {{K,isTorSurjective(R),ideal R}}; --store kSi and surjectivity of Tor maps
		if (not K_0) then print(ideal R);
		); --Appropriate values are now stored
	TorTest := apply(Ts , C -> (
		val := true; --innocent until proven guilty
		for i from 0 to (#(C_0) - 2) do (
			val = val and ((not (C_1)_i) or (not (C_0)_(n+1-i))); --test conjecture
			);
		if not val then (print C#2);
		val --return truth value
		));--rinse and repeat
	{tally(Ks), tally(TorTest)}
	)
--test(ZZ/5,3,4,1)



--For monomial ideals generated at varous degrees
test = method()
test(Ring,ZZ,ZZ,ZZ) := (k,n,mexp,extragens) -> (
	R := ZZ;
	L := {}; --initialize empty list
	I := {}; --ditto
	K := null;
	apply(100, i -> (
		R = randomArtinMonomialRing(k,n,mexp,extragens); --random ring
		K = kSi(R,10); --compute kSi up to some bound
		L = L | {K}; --make a list of the kSi lists
		I = I | {{ideal R, K}}; --make a list of the (I,kSi) pairs
		));
	print(tally(L)); --tally the types
	I
	)


--For monomial ideals generated at varous degrees
--will look at their kSS's and classify
test = method()
test(Ring,ZZ,ZZ,ZZ) := (k,n,mexp,extragens) -> (
	R := null; --arpitrary value
	L := {}; --initialize empty list
	I := {}; --ditto
	K := null;
	apply(100, i -> (
		R = randomArtinMonomialRing(k,n,mexp,extragens); --random ring
		K = kSS(R,10); --compute kSi up to some bound
		L = L | {K}; --make a list of the kSi lists
		I = I | {{ideal R, K}}; --make a list of the (I,kSi) pairs
		));
	print(tally(L)); --tally the types
	I
	)

--time test(ZZ/5,2,5,3)


--For monomial ideals generated at same degree
--checks their kSS
test = method()
test(Ring,ZZ,ZZ,ZZ) := (k,n,d,extragens) -> (
	R := ZZ;
	L := {}; --initialize empty list
	K := null;
	apply(100, i -> (
		R = randomArtinDegreeMonomialRing(k,n,d,extragens); --random ring
		K = kSS(R,10); --compute kSi up to some bound
		L = L | {K}; --make a list of the kSS lists
		));
	tally(L) --tally the types
	)

--time test(ZZ/5,3,2,1)


-*
Takes a ring S, num of gens n, a degree d, and num of
extra generators extragens
*-
test = method()
test(Ring, ZZ, ZZ, ZZ) := (S,n,d,extragens) -> (
	R         := S[vars(0..(n-1))]; --Polynomial Ring
	monos     := monomialList(R,d);
	standard  := {};
	dataList  := {};
	I         := null; --blank ideal
	K         := null; --blank list
	counter   := 0; --counter for debugging purposes
	for i from 0 to (n-1) do (
		standard = standard | {R_i^d};
		);
	combo := toList(0..(extragens - 1));
	while combo#(-1) < #monos do (
		counter = counter + 1; --increment
		--print counter; --display count
		if isReducedList(monos_combo) then (
			print counter;
			I = ideal(standard | monos_combo); --current ideal
			K = kSi(R/I,10); --compute kSS
			dataList = dataList | {K}; --change the second arg in kSS to get loger res
			print I;
			print K;
			);
		combo = nextList(combo);
		);
	tally dataList
	)

time test(ZZ/5,4,3,1)
--try time test(ZZ/5,4,3,2) with the resolution up to length 10 on campus.
--Looks like it will take a while


restart
loadPackage( "SocleSummands", Reload=>true)
    needsPackage "Points";
time test1(ZZ/101,3,3,1)
time test(ZZ/101,3,3,1)
betti F
Rbar**F.dd_3

r=3;n= 6;

F = resOfRandomPoints(3,27);


Rbar = ringOfRandomPoints(2,6)
kSS(Rbar,8)

netList apply(toList(4..9), p->(p,resOfRandomPoints (binomial(p,2)+p//2)))

netList select(apply(toList(4..50), n -> (
	Rbar = ringOfRandomPoints(3,n);
        mat = (res ideal Rbar).dd_3;
    (degree Rbar,betti mat, kSS(Rbar,3)))), t -> t_2 ==t_2)
toList(4..50)

select({1,2,3}, i-> i==3)

n = 12
time I = ideal (Rbar = ringOfRandomPoints(3,n))
betti res I
time kSS(Rbar,7) -- {4,6} ; last matrix of res I has all quadratic terms. Note: missing r+2 is new.

betti res oo
codim I

Rbar
random(1,Rbar)

----
restart
load  "SocleSummands.m2"
debug SocleSummands
S = ZZ/101[a,b,MonomialOrder => {Position => Down}]

m = random(S^{0,1}, S^{-1,-2,-2})
m' = random(S^{0,1}, S^{-2,-2,-2})
I = minors(2,m)
I' = minors(2, m')
R = S/I
R' = S/I'
kSS(R, 8)
kSS2(R,8)
kSS(R', 8)
kSS2(R', 8)
use R'
F = res(coker vars R', LengthLimit=>8)
F.dd_3

--------------
restart
loadPackage "SocleSummands"
debug SocleSummands
installPackage "Points"
S = ZZ/32003[a,b,c,d]
R = S/ideal"a3,b3,c3,d3,cd2"

S = ZZ/32003[a,b,c]
setRandomSeed 0
elapsedTime scan(1, i-> (
	R = S/(ideal apply (gens S, x->x^3) + ideal random(S^1,S^{-3}));
	<<kSS(R,8)<<flush<<endl;
	))
setRandomSeed 0
elapsedTime scan(1, i-> (
	R = S/(ideal apply (gens S, x->x^3) + ideal random(S^1,S^{-3}));
	<<kSS2(R,8)<<flush<<endl;
	))
setRandomSeed 0
elapsedTime scan(1, i-> (
	R = S/(ideal apply (gens S, x->x^3) + ideal random(S^1,S^{-3}));
	<<kSS3(R,8)<<flush<<endl;
	))

time kSS (R,8)
time kSS2 (R,8)
time kSS (R',8)
time kSS2 (R',8)

--Note: these give different answers! something's wrong.
S = ZZ/32003[a,b,c]
R = S/(ideal apply (gens S, x->x^3) + ideal"a2b")
kSS(R,8)
kSS2(R,8)
kSS3(R,8)

use S
R = S/(ideal apply (gens S, x->x^3) + ideal"abc")
kSS(R,8)
kSS2(R,8)
kSS3(R,8)

F = res coker vars R
m = F.dd_3
socleSummand image F.dd_2
kSS (R,10)
kSi (R,8)
socleSummand2 res (coker vars R, LengthLimit=>8)
socleSummand3 res (coker vars R, LengthLimit=>8)
F = res (coker vars R, LengthLimit=>8)
F.dd_2
vars R
trim ideal generalLinearRow F.dd_2
numgens oo

-------
restart
loadPackage "SocleSummands"
debug SocleSummands
S = ZZ/5[a,b,c]
I = ideal"a3,b3,c3"
R = S/I^3
maxres = 7
elapsedTime kSS2(R, maxres)
elapsedTime kSS(R, 8)

----------

///
installPackage "Points"

--The following code leads to the conjecture that the conductor of the 
--homogeneous coordinate ring of n general points in P^r is always a power
--of the maximal ideal. 

--if binomial(r+k-1,r)<n<= binomial(r+k,r) then the conductor is mm^k

--If this is true, then every indecomposable Ulrich module
--over this ring is the factor ring corresponding to a single point.

--Note: for points on the RNC in P^3, the conductor is not a power of mm starting
--with 8 points.

restart
needsPackage "Points"

S = ZZ/32003[a,b,c,d]
R = S/monomialCurveIdeal(S,{1,2,3})
hilbertFunction(2,R)
mm = ideal vars R
I = trim(mm^4+ideal a^5)
isPowerOfMaxIdeal I

pointsOnRNC = method()
pointsOnRNC (Ring,ZZ) := Matrix => (S,n) ->(
    r := numgens S - 1;
    J := apply(n, j-> random(1000));
    M = map(S^(r+1), S^n, (i,j) -> (J_j)^i)
    )

    

isPowerOfMaxIdeal = method()
isPowerOfMaxIdeal Ideal := Boolean => I -> (
    R := ring I;
    degs := I_*/degree;
    e := degs_0;
    (all(degs, d -> d == e) and 
	numgens I == numcols basis(e,ring I)
    )
)

conductorOfPoints = method()
conductorOfPoints Matrix := Ideal => PM -> (
    Rpoints := ring PM/points PM;
    n := numcols PM;
    L := toList(0..n-1);
    partialLists := apply(L, i -> drop(L,-(n-i))|drop(L,i+1));
    C := intersect apply(n, i-> points(PM_{i})+points PM_(partialLists_i));
    trim sub(C,Rpoints)
    )

S = ZZ/32003[a,b,c,d]
R = S/monomialCurveIdeal(S,{1,2,3})
n = 10

scan(30,n->(
	PM = pointsOnRNC(R,n+2);
	if n+2 == degree points PM then
	<<isPowerOfMaxIdeal conductorOfPoints PM)<<endl)

P = pointsOnRNC(R,20)
ring P
degree points P
ring conductorOfPoints P
isPowerOfMaxIdeal conductorOfPoints P

points randomPointsMat(R,4)
degree oo


loadPackage("SocleSummands", Reload=>true)
needsPackage "Points"
S = kk[a,b,c]

scan(30,i-> (R = S/(points pointsOnRNC(S,i+4));
<<degree R<<"  "<<kSS (R,5)<<endl;
))

scan(15,i-> (R = S/(points randomPointsMat(S,i+4));
<<betti res ideal R<<endl;
<<degree R<<"  "<<kSS (R,5)<<endl;)
)

S = kk[a,b,c,d]
scan(21,i-> (R = S/(points pointsOnRNC(S,i+4));
	<<time betti res coker vars R<<endl;
	))


---
restart
needsPackage "SocleSummands"
needsPackage "DGAlgebras"

S = kk[a,b,c]
n = 2
i0 = ideal apply(gens S, x->x^n)
i = i0+ideal"ab"
i = i
m = 2
im = ideal apply(i_*, g->g^m)

kSi (S/im, 8)
isGolod (S/im)

-------------------------
--Golod cases in 3 variables
restart
needsPackage "DGAlgebras"
needsPackage "SocleSummands"
--Long's test for Golod in monomial ideals in 3 vars.
--(1) [I : x1] · [I : (x2, x3)] ⊆ I for all permutations {x1, x2, x3} of {x, y, z}.
--(2) [I : x1] · [I : x2] ⊆ x3[I : (x1, x2)] + I for all permutations {x1, x2, x3} of {x, y, z}.

///
S = ZZ/5[a,b,c]
I = (ideal random(S^1, S^{-2,-2,-1}))*(ideal random(S^1, S^{-2,-1,-1}))
--I = (ideal random(S^1, S^{-2,-2,-1}))
time isGolod (S/I)
time golodTests I

use S
d = 3
F = random(d,S)
G = random(d,S)
I     = ideal(a^2, b^2,a*F+b*G)

G =  ideal(a^2,b^2, a*c,b*c, c^2-a*b)
G =  ideal(a^2,b^2,c^2)
J = apply(II1, I1 -> I1*G);

II1 = apply(10, i-> ideal random(S^1,S^{-2,-3,-4,-5}));
II2 = apply(10, i-> ideal random(S^1,S^{2:-2,2:-4}));
J = flatten apply(II1, I1 -> apply(II2, I2 -> I1*I2));

apply(#J, n -> kSS (S/J_n,5))
apply(#J, n -> minimalBetti J_n)
kSS(S/J_0,5)

scan(#J, n -> <<n<<"  "<<isGolod(S/J_n)<<"  "<<golodTests J_n<<endl)
scan(#II1, n -> <<n<<"  "<<isGolod (S/J_n)<<"  "<<golodTests J_n<<endl)

golodTests(ideal(a^2,b^2, a*c,b*c, c^2-a*b))
use S

kSS(S/G^2,10)
betti res G
R = S/G^2


------------------------------
restart
needsPackage "SocleSummands"
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

S = ZZ/101[a..c]
I0 = monomialIdeal"a3,b3,c4"
--L = orbitRepresentatives(S,I0,{3,3,3,3,3});#L
L = orbitRepresentatives(S,I0,{4,4});#L
LB = select(L, ell -> burchIndex ell >0);#LB
LnotB= select(L, ell -> burchIndex ell ==0);#LnotB
LG = select(L, ell -> isGolod(S/ell));#LG
LGB = select(LG, ell -> burchIndex ell >0);#LGB

L2 = unique apply(subsets(L,2),ell->trim(ell_0*ell_1));#L2

L2B = select(L2, ell -> burchIndex ell >0);#L2B
L2notB= select(L2, ell -> burchIndex ell ==0);#L2notB
L2G = select(L2, ell -> isGolod(S/ell));#L2G
L2GB = select(L2G, ell -> burchIndex ell >0);#L2GB

--

elapsedTime results = apply(L2notB, I ->(
R := S/I;
K:=koszul vars R;
Kcycles := for i from 2 to length K list(
     hasSocleSummand1 image syz K.dd_(i-1)
    );
{isBurch I,Kcycles, isGolod R, kSS(I, 7),I}
));
print"isBurch, KcyclesSummand_2, KcyclesSummand_3, isGolod, kSS, Ideal "
netList results


isGolod(S/I)
kSS(S/I, 9)
kSI(S/I,8)
golodTests I
betti res I

k2 = koszul(2, vars(S/ell))
hasSocleSummand  k2
socle(S^1/ell)
burchIndex ell
kSS(S/ell, 7)
hasSocleSummand coker presentation image k2
socle coker presentation image k2

L = orbitRepresentatives(S,{3,4,4,5});
L1 = unique (L, I-> trim(I + ideal"a3,b3,c3"));
#L
#L1_0

L2 = L1/(I-> S/I)
#L2
time unique(L2/(R->kSS(R,7)))

S = ZZ/101[a..d]
L = orbitRepresentatives(S,{3,3,4});
L1 = unique apply(L, I-> trim(I + ideal"a3,b3,c3,d3"));
#L1
L2 = L1/(I-> S/I);
time tally(L2/(R->kSS(R,7)))

time unique(L2/(R->kSS(R,7)))
time unique(L2/(R->kSS(R,8)))


--generic example
restart
loadPackage "SocleSummands"
n = 5
S = ZZ/101[x_1..x_n]
d = 3
q = sum(n-1,i->S_i*S_(i+1)) -- works with n = 4, d=2; n= 5, d+3--that is, gives kSS == {n-1,n,...}
--q' = sum(n, i->S_i*S_((i+1)%n)) doesn't work with n= 4!
q = S_0*S_1+S_2*S_3 -- works with n=4, d=2
q = S_0*S_1+S_2*S_3 + S_4*S_5 --- n= 6,d=2 gives {3,4,5...}
I = ideal q + ideal apply(n,i->S_i^d)
kSS(S/I,n+2)

kSS(S/I,2*n-2) == {n-1,n,n+1..2*n-3}

K = koszul vars S
d = 2
e = gens K_2*random(K_2, S^{-2-d})
I1 = ideal (K.dd_2*e)
I2 = ideal(e*K.dd_1)
I = trim(I1+I2)
R = S/I
dim R

primaryDecomposition I
kSS(R,3)
sub(e,R)
syz(vars R)
vars R*e

restart
needsPackage "DGAlgebras"
n = 3
kk = ZZ/101
S = kk[x_1..x_n]
d = 3

I = ideal apply(gens S,y->y^d)
I = ideal(random(S^1,S^{n:-d}))
betti res I
m = ideal vars S
e = (n//2)*d+1
J = truncate(e,I)
J' = truncate(e-1,I)
--J' is not Golod (nontrivial homology product); but, conjecture:J is Golod.
betti res J
R = S/J
A = acyclicClosure(R, EndDegree => 0)
isGolod R
isHomologyAlgebraTrivial A

betti res (m*I)
betti res I

restart
needsPackage "DGAlgebras"
kk = ZZ/101
S = kk[a..e]
I = ideal"ab2, cd2, e3,abcd, d2e2, b2e2, ace, b2d2e"
m = ideal vars S
R = S/(m*I)
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/(m*I))
isHomologyAlgebraTrivial A

betti res I
betti res truncate(4,I)
betti res (m*I)

--random examples quadrics and m^3
restart
loadPackage("SocleSummands", Reload => true)
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

isBurch = R -> (
    if not basis(3,R) == 0 then error ("only works for m^3 = 0 rings");
    A := res coker presentation R;
    ty := rank A_4;
    F = res (coker vars R, LengthLimit => 2);
    rank F_2 > (numgens R)^2 - ty
    )

isBurch1 = R -> (
    I = ideal presentation R;
    m = ideal vars ring I;
    0 != (gens (I:m)) % (m*I:m)
	)
    

S = ZZ/3[a..d]
mm = (ideal vars S)
cubes = ideal apply(numgens S, i->  S_i^3)

flatten for m from 5 to 9 list(
for i from 1 to 10 list(
I = trim(ideal random(S^1, S^{m:-2})+cubes);
R = S/I;
ans:= 0;
print isBurch1 R;
if not isGolod R and not isBurch R then(
ans = kSS(R,8);
<<ans<<endl;
ans)))

isGolod R
isBurch R
kSS(R,8)
kSS(R,8)
isBurch R
----------------------
---
restart
needsPackage "DGAlgebras"
n = 3
kk = ZZ/101
S = kk[x_1..x_n]
d = 3
L = {2,4,7}
I = ideal apply(#L,i-> S_i^(L_i))
m = ideal vars S
betti res I
betti res (m*I)

R = S/(m*I)
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/(m*I))
isHomologyAlgebraTrivial A


--Roos example:
restart
needsPackage "DGAlgebras"
kk = ZZ/101
S = kk[a..d]
I = ideal(a^3, a^2*b, (a+d)*(c^2+b^2), c*d^2, d^3)
m = ideal vars S
R = S/(I)
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/(I))
isGolod (S/(m*I))
isHomologyAlgebraTrivial A


---
--Roos example: Claimed to be non-Golod with trivial homology algebra.
restart
needsPackage "DGAlgebras"
needsPackage "SocleSummands"
kk = ZZ/101
--kk = QQ
S = kk[x,y,z,u]
I = ideal(u^3, x*y^2, (x+y)*z^2, x^2*u+z*u^2, y*y*u+x*z*u, y^2*z+y*z^2) -- has the betti nums as in Roos
--I = ideal(u^3, x*y^2, (x+y)*z^2, x^2*u+z*u^2, z^3+x*z*u, y^2*z+y*z^2) -- has the betti nums in

betti (A = res I)
A.dd_1
A.dd_2
R = S/I
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/I)
isHomologyAlgebraTrivial A

betti res( coker vars R, LengthLimit =>7)

T = QQ[t]
((1+t)^4)*sum(10, i-> (6*t^2+12*t^3+9*t^4+2*t^5)^i)

kSS(R,5)
---Roos "exceptional" examples from 
--Homological properties of the homology algebra of the Koszul complex of a local ring:Examples and questions
restart
needsPackage "SocleSummands"
kk = ZZ/101

S = kk[x,y,z,u]
I = ideal"x2, xy, xz + u2, xu, y2 + z2, zu"
betti res I
R = S/I
kSS(R,7) =={0}

use S
I = ideal "xz + u2, xy, xu, x2, y2 +z2,zu,yz"
betti res I
R = S/I
kSS(R, 7) == {0}

use S
I = ideal"x2 + yz + u2,xy,zu,z2,xz + yu,xu"
betti res I
R = S/I
--kSS(R,9) == {5,6,7,8} --this is fairly slow.

-----Roos examples!
needsPackage"DGAlgebras"
kk = ZZ/101
S = kk[u,x,y,z]
use S
II = {I1 = ideal(u^3, x*y^2, (x+y)*z^2, x^2*u+z*u^2, y^2*u+x*z*u, y^2*z+y*z^2),
-- kSS(R,7)={3,4,5}
I2= ideal"u3,y2u2,z2u2,y2zu,xy2,(x+y)z2,x2yz,x2u",
-- kSS(R,7)={0}
I3=ideal"u3,y2u2,z2u2,y2zu,xy2,(x+y)z2,x2yz,x2u+xyu",
-- kSS(R,7)={0}
I4=ideal"u3, y2u2, z2u2, y2zu, xy2, (x + y)z2, x2yz, x2u + zu2",
-- kSS(R,7)={0}
-- the next two has only cubics
I5=ideal"u3, xy2,(x+y)z2, x2u + zu2, y2u + xzu, y2z + yz2",

-- kSS(R,7)={3,4,5}
I6= ideal"u3, x2u, xz2 + yz2, xy2, x2y+y2u, y2z+z2u"
-- kSS(R,7)={4,5,6}
}
netList for I in II list(
H = acyclicClosure(S/I, EndDegree => 0);
isGolod(S/I), isHomologyAlgebraTrivial(H))
elapsedTime kSS(S/I6, 8)

restart
needsPackage "SocleSummands"
needsPackage "DGAlgebras"
kk = ZZ/101
S = kk[x,y,z]
I0 = ideal apply(gens S, x ->x^3)
I = I0+ideal"xyz"
R = S/I
describe S
isGolod(S/I)
describe S
R === S/I
T = torAlgebra(R)
degrees T
describe T
betti res(coker presentation R)
help torAlgebra
betti res I

torAlgebra(S/I0)
--------------------
restart
--needsPackage "StableResidual"
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"
needsPackage "LocalRings"
1,5,5,1
1,5,6,2
1,5,7,3
1,5,8,4
(almost ci, 1,5,5,1)) 
--help MonomialOrbits

generalLink = J -> (
J0 := gens(J)*random(S^(rank source gens J), S^3);
J' := localTrim (ideal(J0):(J))
)

sumBettis = J ->(
    G := localResolution generalLink J;
    sum(4, i-> rank G_i)
    )

dropBettis = J ->(
F := res J;
G := localResolution generalLink J;
sumBettis F> sumBettis G
)

localTrim = J -> ideal (gens J %(M*J))

--Long's test for Golod in monomial ideals in 3 vars; conjecture for all ideals 
--(1) [I : x1] · [I : (x2, x3)] ⊆ I for all permutations {x1, x2, x3} of {x, y, z}.
--(2) [I : x1] · [I : x2] ⊆ x3[I : (x1, x2)] + I for all permutations {x1, x2, x3} of {x, y, z}.

--this shows: no 5 generator Golod monomial ideals (any fin length golod in 3 vars that contains
--x^a,y^b also contains x^(a-1)y^(b-1), and permutations -- two more gens aren't enough.
    

testGolod = I -> (
    P := {{0,1,2},{1,2,0},{2,0,1}};
    G := gens S;    
    t1 := all(3, i-> (
	    g = G_(P_i);
    gens((I:g_0)*(I:ideal(g_1,g_2)))%I == 0
    and
    gens((I:g_0)*(I:g_1))%(g_2*(I:ideal(g_0,g_1))+I) == 0
)))
testGolod ((ideal gens S)^2)    
isGolod(S/I)



S = ZZ/101[x,y,z]
M = ideal vars S
setMaxIdeal ideal vars S

golodExamples = (S, D, E) -> (
I0 := monomialIdeal apply(#gens S, i -> S_i^(D_i));
L := orbitRepresentatives(S, I0, E);
print(#L);
select(L, I-> isGolod (S/I))
)

use S
test = D ->(
I0 = monomialIdeal ideal"x10,y10,z10";
L= orbitRepresentatives(S,I0,D);
--print(#L);
elapsedTime all(L, I -> not testGolod ideal I)) --3 sec

for j from 5 to 10 do <<(j, all(10, i-> test{j,i+2}))<<flush<<endl;


elapsedTime L/(I -> isGolod(S/I)) --23 sec
random 10
golodExamples(S,{10,10,10},{random 10, random 10})
golodExamples(S,{3,3,3},{3,3})

-------
--random examples quadrics and m^3
--installPackage "Posets"
restart
needsPackage "Posets"
loadPackage("SocleSummands", Reload => true)
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

--note: lin7 is defined in SocleSummands.m2 as the list of all (up to permutation) 
--115 linearly presented monomial ideals of degree 7

S = ring lin7_0

addedMonomials = I -> compress((gens I)%ideal (I0 (degree (first (I_*)))_0))
isLinearlyPresented = I -> (
    degs:= last degrees presentation module I;
    1 + (min(I_*/degree))_0 == (max degs)_0
    )
I0 = d-> monomialIdeal trim (monomialIdeal apply(gens S, a -> a^d) + (ideal(x,y))^d+(ideal(x,z))^d+(ideal(y,z))^d)
ex = d -> numcols basis(d,S) - numgens I0 d;
complement Ideal := I -> (
    firstmon = I_0;
    d := (degree firstmon)_0;
    J := basis(d,S);
    monomialIdeal (J % I)
    )
complements = L -> (
    firstmon = L_0_0;
    d := (degree firstmon)_0;
    J := basis (d,S);
    for I in L list monomialIdeal (J % I)
    )
lin7_0
I' =  complement (lin7_0)
sort flatten ((I'_*)/exponents)
lin7Comp = complements lin7;
lin7Comp/numgens

d = 7
ex 7
J = ideal basis (d,S)
elapsedTime L = for i from 0 to ex(d) list orbitRepresentatives(S,I0(d),toList(i:d)); --d=7: 898 sec.

netList (lin = select(flatten for i from 0 to ex list for j from 0 to #L_i-1 list (
    I := L_i_j;
    (isLinearlyPresented I, addedMonomials I)
    ), p -> p_0===true))

lcmLattice monomialIdeal (lin_0_1)
texPoset oo
netList for i from 0 to ex(d) list for j from 0 to #L_i-1 list minimalBetti ideal L_i_j

L7 = L;

elapsedTime lin7 = select(flatten L7, I -> isLinearlyPresented I); --5 SEC
#lin7
elapsedTime lin7Comps = complements lin7;
I = complement lin7_0
value I

restart
needsPackage "Posets"
--needsPackage "Graphics"
--viewHelp Graphics

I = complement lin7_10

point  = ell -> (
    --ell is a list of 3 integers;
    --output is the coordinate point
    d := 1.0*sum ell;
    (d-1.0*ell_0, 1.0*ell_1+.5*ell_2)
    )
    
plotPoints = method()

outerTriangle = d -> (
    --first coord is horizontal offset from lower edge.
    --second coord is vertical offset from center.
    L := flatten apply(d+1, j-> apply(j+1, 
	         i->( -.5*j + 1.0*i,d-j)));
    <<"\\begin{center}"<<endl;
    <<"\\begin{tikzpicture}[scale=1, vertices/.style={draw, fill=red, circle, inner sep=2pt}]"<<endl;
    for i from 0 to #L-1 do(
    <<"\\node [vertices] ("<<i<<") at "<<L_i<<"{$L_i$};"<<endl;
    );
    <<"\\end{tikzpicture}"<<endl;
    <<"\\end{center}"<<endl;
    )
d = 4
outerTriangle d

exponents Ideal := I -> flatten apply(I_*, m -> exponents m)

plotPoints Ideal := I ->plotPoints exponents I

plotPoints List := L->(
    --L is a list of triples {i,j,k} such as that produced by plotPoints Ideal.
<<"\\begin{center}"<<endl;
<<"\\begin{tikzpicture}[scale=1, vertices/.style={draw, fill=black, circle, inner sep=2pt}]"<<endl;
for i from 0 to #L-1 do(
    <<"\\node [vertices] ("<<i<<") at "<<point (L_i)<<"{};"<<endl;
    );
<<"\\end{tikzpicture}"<<endl;
<<"\\end{center}"<<endl;
    )
outerTriangle 7


exponents I

--4-var experiments
restart
needsPackage "Posets"
loadPackage("SocleSummands", Reload => true)
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

kk= ZZ/101
S = kk[x,y,u,v]
d = 8
B = basis(d,S);
expo = {{2,2,2,2},{3,1,2,2},{2,2,3,1}}--lin pres
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I

expo = {{2,2,2,2},{3,1,2,2}}
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I -- lin pres!
--so the fact that the monomials omitted lie in a plane (exp m)_3 ==0
--does not make it "like" the 3 -var case. Not so surprising
--that a "blocker in 3 vars" might not block in 4 vars.

kk= ZZ/101
S = kk[x,y,u]
d = 8
B = basis(d,S);
(e1,e2) = ({2,2,2},{3,1,2})
dist({2,2,2},{3,1,2})
expo = {{2,2,2},{3,1,2}} -- not lin pres
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I

restart
loadPackage("SocleSummands", Reload => true)
kk= ZZ/101
S = kk[x,y,u,v]
d = 8
B = basis(d,S);
expo = {{2,2,2,2},{3,1,2,2},{2,2,3,1}}--lin pres
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I
F = dual res I
betti F

expo = {{2,2,2,2},{3,1,2,2}}
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I -- lin pres!

F = res I
G = res ideal B
e = extend( G,F,map(G_0,F_0,1));
prune coker dual e_4

---
S = ZZ/101[a,b,c,d,e]
I = ideal"ac,ad,bd,be,ce"
betti res I
F = res I
F.dd
(F.dd_2)_{4,3,2,1,0}

restart
loadPackage "SocleSummands"
S = ZZ/101[a,b,c]

for i from 1 to 100 do (
    I := ideal random (S^1, S^{7:-3});
    if isLinearlyPresented I then <<endl<<flush<<minimalBetti ideal random (S^1, S^{7:-3})<<endl;
    )

-------Example of a non-burch ks_3
restart
load "SocleSummands.m2"
S = ZZ/101[x,y]
I =ideal"x4,x2y2,y4"
I =ideal"x8,x4y4,y8"
I = ideal"x3, x2y2, y5"
R=S/I
kSS(R,8)
M = coker random(R^2,R^{4:-2})
MSi(M,8)
F = res coker vars R
hasSocleSummand1 image F.dd_3
b = 15
F = res(coker vars R, LengthLimit=>b)
L = for i from 1 to b list (ns=numcols socleSummands coker F.dd_i;
    <<ns<<endl;
    ns)
for i from 1 to b-1 list L_i - L_(i-1)

kSS(R,6)
k=coker vars R
F = res(k, LengthLimit =>5)
F.dd_4
C = image (F.dd_3)_{0..5}
MSi(C,7)
F.dd_3
F = res(C, LengthLimit =>8)
socleSummands image F.dd_2
socleSummands coker F.dd_3
D = image((F.dd_2)_{0..17})
MSi(D,7)
hasSocleSummand D
for i from 1 to 8 list hasSocleSummand coker F.dd_i
M = image F.dd_2
prune socle M

N = coker random(R^2, R^{4:-2})
elapsedTime MSi(N, 10)
F = res(N, LengthLimit => 8)
elapsedTime for i from 1 to 11 list hasSocleSummand image F.dd_i
use S
mmS = ideal gens S
I = ideal"x2,y2"
J =trim( mmS*I + ideal"x2, xy")
J = ideal"x2,xy,y4"
(J:mmS)/J
presentation(I/J)
B = S/J
A = S/I
AB = map(A,B)
presentation (module ideal vars A)
M = image oo
socleSummands module ideal gens B
M' = image (vars B)_{1}
MSi(M', 6)
F = res M'
F.dd_1
M'' = image (vars A)_{1}
MSi(M'', 6)
F = res(coker vars B, LengthLimit=>10)
L = for i from 1 to 10 list numcols socleSummands coker F.dd_i
for i from 1 to 9 list L_i - L_(i-1)
 ----------
 
--What is kSS for n random points in P^2 when they are not burch?
restart
loadPackage "SocleSummands"
loadPackage "Points"
p = 7
n = binomial(p,2) + p//2
ringOfRandomPoints(11,2)

test = n->(
IS = points randomPointsMat(S,n);
I = TS IS;
<< n << "  " <<kSS(T/I, 7)<<endl)
for n from 4 to 40 do test n





test n
n=5
kk = ZZ/32003
--kk = GF(32)

restart
load "SocleSummands.m2"
kk = ZZ/101
S = kk[a, x,y,z]
I = ideal(a^3*x, y^4,z^4, a^3*y - z^3*y, a^3*z - y^3*x)
R = S/I
mmR = ideal gens R
J = ideal(x+y)+mmR^2

MSi(module J,8)


needsPackage SocleSummands
S = kk[x,y,z]
kSS (S/I, 9)
res I
mm
T = kk[x,y]
I = ideal(x^3, x^2*y^2, y^3) 
I = ideal(x^4, x^2*y^2, y^4)
R = T/I
mmR = ideal gens R
J = ideal(x+y)+mmR^2
MSi(module J, 8)

TS = map(T,S,random(T^1,T^{3:-1}));
matrix TS




needsPackage "Points"
S
n=11
IS = points randomPointsMat(S,n);
R = S/IS
mmR^3 = ideal gens R

weakLin(module (mmR^2), 7)

 n-1 then << toString I <<endl<<endl;
    )

n = 3
S = ZZ/101[x_0..x_(n-1)]
T = ZZ/101[x_0..x_1]
TS = map(T,S,random(T^1, T^{n:-1}))
mm = monomialIdeal gens S
d = 5
ze = monomialIdeal(0_S)
L = orbitRepresentatives(S,ideal ze, mm^10, 6);
elapsedTime L = orbitRepresentatives(S,toList(4:4)|toList(7:5));

#L

I = (mm^[2])^2
I = L_2
J = (I*mm):(I:mm);
numcols compress ((gens J) % mm^2)
elapsedTime (L/test)

I = monomialIdeal(x_0^4,x_0^3*x_1,x_0^2*x_1^2,x_1^5,x_1^3*x_2,x_0^3*x_2^2,x_0*x_1^2*x_2^2,x_0^2*x_2^3,x_0*x_1*x_2^3,x_1*x_2^4,x_2^5)
R = S/I
use R
S

data I    

I = ideal"x2,xy,y2"    
data I

I = ideal"x2, xy2, y3"
data I

S = kk[a,b]
I=ideal(a^4, a^2*b, b^2)
data I


loadPackage "SocleSummands"
Lburch3

I = Lburch1_5
S = ring I
JS = (I*mm) : (I:mm)
R = S/I
J = sub(JS,R)
weakLin(J, 7)

MSi(module J,5), 


installPackage "CorrespondenceScrolls"
viewHelp carpet
I = carpet{2,7}
betti res I
hilbertSeries((ring I)^1/I, Reduce =>true)
S = ring I
vars S
kk = coefficientRing S
T = kk[y_0..y_6]
TS = map(T,S,random(T^1, T^{11:-1}))
J = TS I
hilbertSeries(T^1/J, Reduce =>true)

restart
kk = ZZ/101
m=3
n=5
d = 2
S = kk[x_(1,1)..x_(m,n)]
M = genericMatrix(S,x_(1,1),m,n)
J = minors(d,M)
T = kk[y_1..y_((m-d+1)*(n-d+1))]
TS = map(T,S,random(T^1,T^{m*n:-1}))
I  = TS J
mm = ideal gens T
(mm*I) :(I:mm)

------
restart
loadPackage "SocleSummands"
needsPackage "Points"
kk = ZZ/101
S = kk[x,y,z]
I = points randomPointsMat(S, 5)
M = syz gens I
R = S/I
N = sub(M,R)
J = coker transpose (N_{0})
J = coker random(R^1, R^{2:-2})
J = coker matrix"x2,y2"
weakLin(J, 7)
weakLin(J, 7,2)
trim (ideal gens R)^2



F = res J
sy = F.dd_3
sy = coker F.dd_4
prune(H =  Hom(sy,sy))
ho = apply (numgens H, i ->homomorphism H_{i})
T = kk[a_0..a_3]
hoT = apply(ho, A -> sub(matrix A, T))
HoT = sum(4, i-> hoT_i*a_i)
codim minors(2,HoT)
degree minors(2, HoT)
F.dd

vars S

restart
loadPackage "SocleSummands"
needsPackage "AInfinity"
kk = ZZ/101
U = kk[a,b]
m = ideal gens U
I = m^[3]*a*b
I = (ideal"a2,b2")^2
I = (ideal"a3,b3")^2
R = U/I
res coker vars R
mR = sub(m, R)
weakLin(coker random(R^1, R^{-3}),7)
--!
betti res(coker random(R^1, R^{-3}), LengthLimit => 7)
weakLin(module ideal "a3+b3",7)
MSi(module ideal "a3+b3", 7)
MSi(module (mR^2),7)
weakLin(module (mR^2),5)
betti burkeResolution(module ideal "a3+b3", 7)
betti (F = burkeResolution(coker random(R^1, R^{-3}), 7))
picture F
use U
betti res I
RU = map(R,U)
betti (G = res (pushForward(RU, coker random(R^1, R^{-3}))))
G.dd_2

------
needsPackage "MonomialOrbits"
needsPackage "SocleSummands"
kk = ZZ/101
S = kk[a,b,c]
m = ideal gens S
L = orbitRepresentatives(S, m^[5], {4})
netList apply(L, I -> weakLin (S^1/I,5))
I = last L
weakLin(S^1/(last L), 5)
ids = apply(L,I -> (
	F := res I;
	{ideal F.dd_1, ideal F.dd_2, ideal F.dd_3}
	));
for i from 0 to #ids-1 do(
    (gens iid_i_0+ids_i_1  

R = S/m^[3]
J = ideal(a^2*b^2)
weakLin(module J, 7)
m*L_0 : (L_0:m)

-----------
restart
needsPackage "MonomialOrbits"
needsPackage "SocleSummands"
kk = ZZ/101
S = kk[a..f]
m = ideal gens S
L = orbitRepresentatives(S, m^[5], {4})
I = (m^[3])^2
R = S/I
f = a+b+c+d
weakLin (module ideal (f^2),3)
--!
S = kk[x_(0,0)..x_(2,3)]
A = genericMatrix(S, x_(0,0),3,4)
I = minors(3, A)
R = S/I
N = coker random(R^1,R^{-1,-1,-2})
weakLin(N, 4)

restart
loadPackage ("SocleSummands", Reload =>true)
needsPackage "MonomialOrbits"
kk = ZZ/101
c = 2
d = 2

S = kk[x_0..x_(c-1)]
m = ideal gens S
squareFree = (S,d) ->(
    y := symbol y;
    T = kk[y_0..y_(numgens S -1), SkewCommutative => true];
    (map(S,T,vars S))(ideal basis (d,T))
    )
n = squareFree(S,d)
orbitRepresentatives(S,m^[d],m^d,1)


R = S/m^5
F = res coker vars R
prune socle image oo.dd_2
socleSummands coker F.dd_3


--canonical blowup -- when Gorenstein?
installPackage "ReesAlgebra"
installPackage "Points"

restart
needsPackage "ReesAlgebra"
needsPackage "Points"
needsPackage "FastMinors"
--viewHelp

--viewHelp omegaPoints
--viewHelp Points
canonicalIdeal = method()
canonicalIdeal Ideal := Ideal => I ->(
    S := ring I;
    R := S/I;
    F := res I;
    omega := coker sub(transpose F.dd_(length F), R);
    H := Hom(omega,R^1);
    deg := max degrees source gens H;
    g := (gens H)*random(source gens H, R^-deg);
    trim sub(ideal g,R) ---canonical ideal of a 1-dim ring.
)

kk=ZZ/101
S = kk[x,y,z]

I = points randomPointsMat(S,6)
I = minors(2, random(S^3, S^{-1,-2}))


degree I
betti gens I
R = S/I
om = canonicalIdeal I
codim om
RI = reesIdeal om

betti gens om
T = kk[x,y,z, vars(0..numgens ring RI -1)]
vars ring RI
omegaRI = ker map(ring RI/RI, T, vars ring RI|vars R)
betti (F =res omegaRI)
length F
1+numgens ring RI == length F
codim ideal F.dd_(length F)
E = Ext^2(omegaRI,ring omegaRI)
betti prune E
dim T
codim minors(2, presentation E)
--not Gorenstein?

codim chooseGoodMinors(1000, 2, presentation E)
betti prune E
betti res prune E

---weakLin of modules over 0-dim Gorenstein.
restart
loadPackage "SocleSummands"
kk = ZZ/30223
S = kk[x_1..x_4]
I = trim ideal fromDual(sum(4, i-> S_i^2))
betti (F = res I)

R = S/I
f = ideal random(R^1, R^{-2})
betti (F = res f)
ent1 F
weakLin((ring f)^1/f, 5)
M = coker transpose F.dd_4
betti (G = res M)
ent1 G
--!

restart
loadPackage "SocleSummands"
kk = ZZ/30223
S = kk[x_1..x_3]
I = trim ideal fromDual(sum(numgens S, i-> S_i^3))
betti (F = res I)

R = S/I
RS = map(R,S)
M = coker random(R^1, R^{-2})
betti (F = res M)
dim M
ent1 F
weakLin(M, 5)

needsPackage "AInfinity"
--viewHelp picture
mR = aInfinity R;
m =  aInfinity M;

G = burkeResolution(M, 6)
netList for i from 3 to 6 list picture G.dd_i
netList (Km = sort select(keys (m = aInfinity(M)) , k -> class k_0 === ZZ and sum k <5))
displayBlocks G.dd_5
betti res I

--- 2 vars
restart
loadPackage "SocleSummands"
kk = ZZ/30223
S = kk[x_1..x_2]
I = trim ideal fromDual(sum(numgens S, i-> S_i^3))
betti (F = res I)

R = S/I
RS = map(R,S)
M0 = coker (R_0^2)
betti (F = res (M0, LengthLimit => 6))
MM = for i from 1 to 5 list coker transpose F.dd_(i);

FF = for i from 0 to 4 list res(MM_i, LengthLimit => 5)
netList for i from 0 to 4 list (i, betti FF_i, ent1 FF_i)
i = 2
(i, betti FF_i, ent1 FF_i)
G2 = burkeResolution (MM_2, 5)
displayBlocks G2.dd_3
picture G2.dd_3

N = coker transpose F.dd_2
Gmin = res (N, LengthLimit =>5)
netList for i from 1 to length Gmin list  ent1 G.dd_i
netList for i from 1 to length Gmin list  G.dd_i
needsPackage "AInfinity"
--viewHelp picture
mR = aInfinity R;
m =  aInfinity N;

H = burkeResolution(N, 6)
netList for i from 3 to 6 list picture H.dd_i
netList (Km = sort select(keys (m = aInfinity(N)) , k -> class k_0 === ZZ and sum k < 7))
apply(Km, k->(k,ent1 m#k))

netList for i from 1 to 4 list displayBlocks H.dd_i
betti res I

-----
--1-dimensional minimal mult rings

restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z] -- Degrees => {3,7,8}]
T = kk[t]
I = ker map(T,S, {t^3, t^7, t^8})
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.
mm = ideal vars R
F = res mm
M = coker syz vars R
M' = coker syz matrix"x2,y,z"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2
(x^2*M) : M
(x^2*M') : M'
isIsomorphic(M'1, M'++M', Homogeneous => false) --gives error in quasi-hom case.
isIsomorphic(M1,M++M, Homogeneous => false)
(isIsomorphic(M,M1,Homogeneous =>false))_0 == false

-------
restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z,Degrees => {3,8,10}]
S = kk[x,y,z] --,Degrees => {3,8,10}]
T = kk[t]
I = ker map(T,S, {t^3, t^8, t^10})
betti res (S^1/I)
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.

M = coker syz vars R
M' = coker syz matrix"x2,y,z"
M'' = coker syz matrix"x3,y,z"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2
M''1 = coker (res M'').dd_2

isIsomorphic(M'1, M'++M', Homogeneous => false) --gives error in quasi-hom case.
isIsomorphic(M1,M++M, Homogeneous => false)
isIsomorphic(M''1,M'++M'', Homogeneous => false)
(isIsomorphic(M,M',Homogeneous =>false))_0 == false

-------

restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z]
T = kk[t]
I = ker map(T,S, {t^3, t^10, t^11})
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.

M = coker syz vars R
M' = coker syz matrix"x2,y,z"
M'' = coker syz matrix"x3,y,z"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2
M''1 = coker (res M'').dd_2

isIsomorphic(M1,M++M, Homogeneous => false)
isIsomorphic(M'1, M'++M', Homogeneous => false) --gives error in quasi-hom case.
isIsomorphic(M''1,M''++M'', Homogeneous => false)
(isIsomorphic(M,M',Homogeneous =>false))_0 == false

-------
--arf: m
restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z,u]
T = kk[t]
I = ker map(T,S, {t^4, t^7, t^9, t^10})
--not arf
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.

M = coker syz vars R
M' = coker syz matrix"x2,y,z,u"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2

isIsomorphic(M1,M++M++M, Homogeneous => false)
isIsomorphic(M'1, M++M'++M', Homogeneous => false)
isIsomorphic(M''1,M''++M'', Homogeneous => false)
(isIsomorphic(M,M',Homogeneous =>false))_0 == false

-------
restart
loadPackage "SocleSummands"
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"

kk = ZZ/32003
S = kk[a,b,c]

i = ideal "a4,b4,c4, abc"
kSS(S/i, 8) 
i = ideal "a4,b4,c4, ab3, b2c2"
kSS(S/i, 7) 
i = ideal "a4,b4,c4, a2b2, b2c2"
elapsedTime kSS(S/i, 10) 
i = ideal "a4,b4,c4, ab3,bc3"
isGolod(S/i)
use S
elapsedTime kSS(S/i, 8) 
i = ideal "a4,b4,c4, ab3,b3c"
elapsedTime kSS(S/i, 8) 
i = ideal "a4, a2b2, a2c2"
elapsedTime kSS(S/i, 8) 
isGolod (S/i)

use S
i = ideal"a2,b2"*ideal "a2,b2,c2"
elapsedTime kSS(S/i, 8) 
isGolod (S/i)
use S
golodTests i

use S
I = monomialIdeal "a6,b6,c6"
elapsedTime L = orbitRepresentatives(S, I, {5,5,5,5});#L
LG = select(L, J -> golodTests J and 
                    burchIndex J == 0 );#LG
elapsedTime LG456 = select(LG, J -> kSS(S/J, 6) == {4,5});#LG456		
elapsedTime LG345 = select(LG, J -> kSS(S/J, 4) == {3});#LG345

scan(LG345, I ->(
	R = S/I;
	<<hasSocleSummand ker koszul(2, vars R)<<endl
	)
)
I = LG345_0

LGlin = select(LG, J -> linearFormsInRes J =={true,false})
		    ;#LG
for I in LG do <<I<<" "<< kSS(S/I, 6)<< linearFormsInRes I<<endl
elapsedTime LGlist = for I in LG list{I,kSS(S/I,6), linearFormsInRes I}
select(LGlist, ell -> ell_1!={3,4,5})
LGlist_0
kSS(S/LG_55,8)						    -- 

I = ideal"a3,b3,c3"*ideal"a2,b2,c2"
trim I
burchIndex I
isGolod(S/I)
linearFormsInRes I
kSS(S/I,7)

use S
I = ideal"a3,b3,c3,abc"
trim I
burchIndex I
isGolod(S/I)
linearFormsInRes I
kSS(S/I,7)

kSS(S/LG_0, 8)
F = res I
F.dd_3
golodTests I
F = res LG_0
F.dd_2
R = S/LG_0
G = res coker vars R
G.dd_2

viewHelp EagonResolution
use S
I = ideal"a6,b6,c6, a3b2, b3c2, c3a2"
linearFormsInRes I
kSS(S/I, 7) == {3,4,5}
needsPackage "EagonResolution"
R= S/I

F = eagonResolution(R,5)
picture(F.dd_3,Display => "DisplayBlocks")
mapComponent(F.dd_4,(3,{}),(0,{3}))
picture(F, Display => DisplayBlocks)


---------------
restart
S = ZZ/32003[a_0..a_5]
M = genericMatrix(S, a_0, 2,3)
I = minors(2, M)
R = S/I
M = coker matrix{apply(6, i-> R_i^2)}
F = res(M, LengthLimit => 6);
betti F
tally((flatten entries F.dd_6)/degree)
M1 = submatrixByDegrees(F.dd_5,(7,7),(8,8));

M1 = submatrixByDegrees(F.dd_4,6,7);
betti (G = res (coker M1))
N = coker G.dd_2
M2 = matrix{apply (3, i-> R_i)}
G2 = res coker M2
N2 = coker G2.dd_5
isIsomorphic(N2, N)
betti res N
betti res N2
betti res coker random(R^{1,0}, R^{6:-2})

S = ZZ/32003[a,b]
R = S/(ideal vars S)^3
M = random(R^1,R^{-2})
betti res(coker M,LengthLimit => 10)

-----------
needsPackage "DGAlgebras"
loadPackage"SocleSummands"
netList (L = for i from 1 to 30 list (
    I = ideal ringOfRandomPoints(3,i);
    if burchIndex I == 0 then I else continue
    ))
for I in L do (<<kSS((ring I)/I, 6)<<endl<<flush)
betti res L_2
isGolod (ring L_2/L_2)
degree(ring L_2/L_2)
netList (L = for i from 1 to 20 list (
    I = ideal ringOfRandomPoints(3,i);
    kSS(ring I/I, 6)))
R = ringOfRandomPoints(3,13)
kSS(R,6)
use R


----- in 3 vars: does last map linear imply 
--the second Koszul cycles have socleSummand?
restart
loadPackage "SocleSummands"
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"

kk = ZZ/32003
S = kk[a,b,c]
I = monomialIdeal "a3,b3,c3"
elapsedTime L = orbitRepresentatives(S, I, {5,5,5});#L
LG = select(L, J -> (not golodTests J) and 
                    burchIndex J > 1) ;#LG
elapsedTime LG456 = select(LG, J -> kSS(S/J, 6) == {4,5});#LG456		
elapsedTime LG345 = select(LG, J -> kSS(S/J, 4) == {3});#LG345

J = LG_0
B = drop(flatten entries basis R, 1)
Blin = for f in B list     needsPackage "Points";
    x:= symbol x;
    S := ZZ/32003[x_0..x_r];
    I := points randomPointsMat(S,n);
    R := coefficientRing S[(gens S)_{1..r}];
    phi := map(R,S,{random(1,R), R_0..R_(r-1)});
    assert(codim phi I == r);
    R/phi I)


resOfRandomPoints = (r,n) -> (
    Rbar := ringOfRandomPoints(r,n);
    S = ambient Rbar;
    use S;
    << minimalBetti (S^1/ideal Rbar)<<", "<<betti (F := res(coker vars Rbar, LengthLimit => r+1))<<endl;
    F
    )
linearFormsInRes = I ->(
    mm := ideal vars ring I;
    F := res module I;
    apply(length F , i-> t :=(0 != compress ((gens ideal F.dd_(i+1))%mm^2)))
)

linearDimInRes = I ->(
   mm := ideal vars ring I;
   F := res module I;
   apply(length F , i-> 
       numcols compress (mingens ideal F.dd_(i+1))%mm^2))
      

beginDocumentation()
--Documentation
-*
restart
loadPackage("SocleSummands", Reload=> true)
installPackage "SocleSummands"
check SocleSummands
Todo list:
	"randomMonomial",
	"randomArtinIdeal",
	"randomArtinDegreeIdeal",
	"randomArtinMonomialIdeal",
	"randomArtinDegreeMonomialIdeal",
	"randomArtinRing",
	"randomArtinDegreeRing",
	"randomArtinMonomialRing",
	"randomArtinDegreeMonomialRing",
*-
doc ///
	Key
		hasSocleSummand
		(hasSocleSummand, Module)
	Headline
		Determines if residue field is a direct summand of input
	Usage
		 hasSocleSummand M
	Inputs
		M:Module
			A module over an Artin local ring
	Outputs
		:Boolean
			Whether the residue field is a direct summand of M
	Description
		Text
			Takes in a module over a local Artin ring and determines
			is the residule field of the ring is a direct summand of
			M
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			M = coker matrix {{a*b,b*c},{c,a*b*c}}
			hasSocleSummand M

///

doc ///
	Key
		randomArtinIdeal
		(randomArtinIdeal, Ring, ZZ, ZZ)
		(randomArtinIdeal, Ring, ZZ)
	Headline
		Produces a random Artinian ideal
	Usage
		randomArtinIdeal(R, m, n)
		randomArtinIdeal(R, n)
	Inputs
		R:Ring
		m:ZZ
			that bounds the degree of the generators
		n:ZZ
			that controls the number of minimal generators
	Outputs
		:Ideal
		    with a minimal set of generators
		    
	Description
		Text
			Let r be the number of generators of R. The function 
			constructs a random homogenous Artinian ideal in R with 
			at most r + n generators. The degree of each generator 
			is at most m + 1. If no value of m is given, the 
			default value is 20.
			
			Note that the generators of the ideal includes random powers of 
			the generators of R.
		Example
			S = ZZ/32003[x,y,z,w]
			I = randomArtinIdeal(S,6,4)
			numgens I
			L = flatten entries mingens I;
			apply(L, v -> first degree v)
    	SeeAlso
	     randomMonomial
	     randomArtinDegreeIdeal
	     randomArtinMonomialIdeal
	     randomArtinDegreeMonomialIdeal
	     randomArtinRing
	     randomArtinDegreeRing
	     randomArtinMonomialRing
	     randomArtinDegreeMonomialRing
    
///

doc ///
	Key
		kSi
		(kSi, QuotientRing)
		(kSi, QuotientRing, ZZ)
	Headline
		A boolean list of when k is a direct summand of its syzygy
	Usage
		kSi R
		kSi(R, n)
	Inputs
		R:QuotientRing
			An Atrin local ring
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:List
			List of booleans when k is a direct summand of its
			syzygy
	Description
		Text
			Takes in an Artin local ring, R, and computes a free
			resolution of its residue field, k, and determines when
			k is a direct summand of its ith syzygy. Outputs a list
			of boolean values for when it is.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			elapsedTime kSi(R, 8)
///


doc ///
	Key
		kSS
		(kSS, QuotientRing, ZZ)
	Headline
		A minimial list of generators for kSS(R)
	Usage
		kSS(R, n)
	Inputs
		R:QuotientRing
			An Atrin local ring
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:List
			List of generators for the numerical semigroup generated
			by the indexes for which k is a direct summand of its
			syzygies of the corresponding index
	Description
		Text
			Takes in an Artin local ring, R, and computes a free
			resolution of its residue field, k, and determines when
			k is a direct summand of its ith syzygy. For the idexes
			that return true we return the sublist of minimal
			generators for the numerical semigroup they generate
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			kSS(R, 8)
///


doc ///
	Key
		kSI
		(kSI, QuotientRing, ZZ)
	Headline
		The (kS) index of R
	Usage
		kSI(R, n)
	Inputs
		R:QuotientRing
			An Atrin local ring
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:ZZ
			the infimum of indexes i for which R is (kSi) for i > 0
	Description
		Text
			Takes in an Artin local ring, R, and integer n and
			computes a free resolution of its residue field, k, up
			to index n and determines when k is a direct summand of
			its ith syzygy. If any positive index results in (kSi)
			being true then this will return the minimum of these
			indexes. Otherwise this method returns infinity.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3,a2c2")
			R = S/I
			kSI(R, 8)
///


doc ///
	Key
		MSi
		(MSi, Module, ZZ)
	Headline
		A boolean list of when k is a direct summand of the syzygies of M
	Usage
		MSi(M, n)
	Inputs
		M:Module
			A module over an Atrin local ring, R
		n:ZZ
			maximum index of syzygy module considered
	Outputs
		:List
			List of booleans when k is a direct summand of the
			syzygies of M
	Description
		Text
			Takes in a module over an Artin local ring, R, and
			computes a free resolution of its residue field, k, and
			determines when k is a direct summand of its ith
			syzygy. Outputs a list of boolean values for when it is.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			M = coker matrix {{a*b,b*c},{c,a*b*c}}
			MSi(M, 8)
///

doc ///
	Key
		TorMap
		(TorMap, QuotientRing, ZZ)
		(TorMap, Matrix, Module, ZZ)
	Headline
		Computes the induced map on Tor_i
	Usage
		TorMap(R, i)
		TorMap(f, P, i)
	Inputs
		R:QuotientRing
			A module over an Atrin local ring, R
		f:Matrix
			Module map from M to N
		P:Module
			Fuctor applied is Tor_i(-,P)
		i:ZZ
			Index of Tor
	Outputs
		:Matrix
			The induced map Tor_i(M,P) -> Tor_i(N,P)
	Description
		Text
			Given an Artin local ring R with ambient ring S
			and defining ideal I in S and integer i, TorMap(R,i)
			computes the induced map Tor_i(mI,k) -> Tor_i(I,k)
			where m is the maximal ideal generated by the variables
			of S and k is the residue field of R thought of as an
			S-module
			
			Given modules M,N,P over an Artin local ring R, an
			integer i, and a matrix f:M->N thought of as a module
			map, TorMap(f,P,i) computes the induced map
			Tor_i(M,P) -> Tor_i(N,P)
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			M = R^3
			N = coker matrix {{a*b},{b*c}}
			P = coker matrix {{a},{b},{c}}
			F = matrix {{a^2,b^2,c^2},{a*b,b*c,c*a}}
			f = inducedMap(N,M,F)
			TorMap(f,P,3)
			TorMap(R,3)
///

doc ///
	Key
		TorMaps
		(TorMaps, Matrix, Module, ZZ)
	Headline
		Computes a list of the induced maps on Tor_i
	Usage
		TorMaps(f, P, iend)
	Inputs
		f:Matrix
			Module map from M to N
		P:Module
			Fuctor applied is Tor_i(-,P)
		iend:ZZ
			Maximal index of Tor
	Outputs
		:List
			List of the induced maps Tor_i(M,P) -> Tor_i(N,P) for
			i from 0 to iend
	Description
		Text
			Given modules M,N,P over a ring R, an integer i, and a
			module map f:M->N, TorMaps(f,P,i) computes a list of
			the induced maps Tor_i(M,P) -> Tor_i(N,P) for i from 0 
			to iend
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			M = R^3
			N = coker matrix {{a*b},{b*c}}
			P = coker matrix {{a},{b},{c}}
			F = matrix {{a^2,b^2,c^2},{a*b,b*c,c*a}}
			f = inducedMap(N,M,F)
			TorMaps(f,P,3)
///


doc ///
	Key
		isTorSurjective
		(isTorSurjective, QuotientRing)
		(isTorSurjective, QuotientRing, ZZ)
	Headline
		List surjectivity status of the induce Tor maps
	Usage
		isTorSurjective(R)
		isTorSurjective(R, n)
	Inputs
		R:QuotientRing
			Specified ring
		n:ZZ
			Maximal index for which surjectivity status will be
			computed
	Outputs
		:List
			List of the surjectivity status of the induced maps 
			Tor_i(M,P) -> Tor_i(N,P) for i from 0 to n
	Description
		Text
			For a local Artin ring R with ambient polynomial ring S
			and defining ideal I in S. Let m be the maximal ideal of
			S generated by the variables, the preimage of the
			maximal ideal of R. As S-modules we have the inclusion
			map mI -> I which induces maps Tor_i(M,P) -> Tor_i(N,P).
			This method computes a list booleans for i from 0 to n
			determining if the induced map is surjective.
		Example
			S = ZZ/1993[a,b,c]
			I = ideal("a4,b4,c4,ab3")
			R = S/I
			isTorSurjective(R,3)
///

test1 = method()
test1(Ring, ZZ, ZZ, ZZ) := (k,n,d,extragens) -> (
	S         := k[vars(0..(n-1))]; --Polynomial Ring
	monos     := monomialList(S,d);
	standard  := {};
	dataList  := {};
	I         := null; --blank ideal
	J         := null; --blank ideal	
	K         := null; --blank list
	R         := null;
	Rbar      := null;
	counter   := 0; --counter for debugging purposes
	for i from 0 to (n-1) do (
		standard = standard | {S_i^d};
		);
	combo := toList(0..(extragens - 1));
	while combo#(-1) < #monos do (
		counter = counter + 1; --increment
		--print counter; --display count
		if isReducedList(monos_combo)and counter==1 then (
			print counter;
			I = ideal(standard | monos_combo); --current ideal
			J = I:ideal vars S;
			R = S/I;
			Rbar = R/sub(J,R);
			F := res(coker vars R, LengthLimit => 10);
			--print("okay1");
			
			K = hasSocleSummand1(Rbar, F); --compute kSS
			--print("okay2");
			dataList = dataList | {K}; --change the second arg in kSS to get loger res
			print I;
			print K;
			);
		combo = nextList(combo);
		);
	tally dataList
	)

golodTests = method(Options =>{Verbose => false})
golodTests Ideal :=  o ->I ->  (
    --by paper of Dao and diStefani, this pair of tests is necessary for Golodness,
    --and sufficient for monomial ideals in 3 vars; not known whether its
    --its sufficient for general homogeneous ideals in 3 vars.
    R := ring I;
    v := gens R;
    V := (permutations v)_{0,3,4}; -- the cyclic permutations.
    condition1 := for w in V list(
	gens((I:w_0)*(I:(ideal(w_1,w_2))))%I == 0
	);

    condition2 := for w in V list(
	gens((I:w_0)*(I:w_1)) % (w_2*(I:(ideal(w_0,w_1)))+I) == 0
	);
    a1 := all(condition1, t->t);
    if a1 == false then
        if o.Verbose == false then return false else
	     <<netList I_*<<endl<<"failed condition 1"<<endl;
    a2 := all(condition2, t->t);
    if a2 == false then  
        if o.Verbose == false then return false else 
	      <<netList I_*<<endl<<"failed condition 2"<<endl;    
    true
)    

///
restart
loadPackage("SocleSummands", Reload => true)
R = ZZ/101[a,b,c]
I = ideal"a2,b3,c4,ab,bc"
golodTests I

restart
uninstallPackage "SocleSummands"
time installPackage "SocleSummands"
--Checko "Soclesummands"
needsPackage "SocleSummands"
viewHelp SocleSummands
time test1(ZZ/101, 3,3,1)
///							    -- 

///
dist1 = (e1,e2) -> sum apply(#e1, i-> max (e1_i,e2_i)) - sum e1
--e1,e2 should be lists of ZZ of the same length and same sum d, representing two monomials of degree d.
--the distance is the degree of the lcm minus d.
-- this is not exported.

dist = method()
dist(List, List) := ZZ => (E1,E2) -> (
    --E1,E2 should each be a list of lists of ZZ, representing monomials of the same degree.
    --the distance is the minimal distance between monomials in the set.
    min flatten apply(E1, e1-> apply(E2, e2-> dist1(e1,e2)))
	)
///

-*
returns a list of n complexes of length ell:
K_0 = koszul vars R
K_1 = result of killing cycles of degree 1
K_2 = result of killing cycles of degree 1, 2
etc
*-

restart
loadPackage("SocleSummands", Reload=> true)
uninstallPackage "SocleSummands"
restart
installPackage "SocleSummands"
--check "SocleSummands"




--test on monomial ideals gen at same degree
--randomArtinDegreeMonomialRing(R,n,d,extragens)
--Things to test
-*
1) (kSi) on the rings
2) The Tor^S_i surjectivity for n+1-i > 0, i < n+1

--Tally the (kSi) patterns
--verify tor surjectivity conjecture
*-
test = method()
test(Ring, ZZ, ZZ, ZZ) := (k,n,d,extragens) -> (
	Ks := {}; --Store (kSi)
	K  := null;
	Ts := {}; --Store surjectiveity of ith Tor map
	R  := ZZ;
	for i from 1 to 100 do (
		R = randomArtinDegreeMonomialRing(k,n,d,extragens); --random ring
		K  = kSi(R);
		Ks = Ks | {K};  --store kSi
		Ts = Ts | {{K,isTorSurjective(R),ideal R}}; --store kSi and surjectivity of Tor maps
		if (not K_0) then print(ideal R);
		); --Appropriate values are now stored
	TorTest := apply(Ts , C -> (
		val := true; --innocent until proven guilty
		for i from 0 to (#(C_0) - 2) do (
			val = val and ((not (C_1)_i) or (not (C_0)_(n+1-i))); --test conjecture
			);
		if not val then (print C#2);
		val --return truth value
		));--rinse and repeat
	{tally(Ks), tally(TorTest)}
	)
--test(ZZ/5,3,4,1)



--For monomial ideals generated at varous degrees
test = method()
test(Ring,ZZ,ZZ,ZZ) := (k,n,mexp,extragens) -> (
	R := ZZ;
	L := {}; --initialize empty list
	I := {}; --ditto
	K := null;
	apply(100, i -> (
		R = randomArtinMonomialRing(k,n,mexp,extragens); --random ring
		K = kSi(R,10); --compute kSi up to some bound
		L = L | {K}; --make a list of the kSi lists
		I = I | {{ideal R, K}}; --make a list of the (I,kSi) pairs
		));
	print(tally(L)); --tally the types
	I
	)


--For monomial ideals generated at varous degrees
--will look at their kSS's and classify
test = method()
test(Ring,ZZ,ZZ,ZZ) := (k,n,mexp,extragens) -> (
	R := null; --arpitrary value
	L := {}; --initialize empty list
	I := {}; --ditto
	K := null;
	apply(100, i -> (
		R = randomArtinMonomialRing(k,n,mexp,extragens); --random ring
		K = kSS(R,10); --compute kSi up to some bound
		L = L | {K}; --make a list of the kSi lists
		I = I | {{ideal R, K}}; --make a list of the (I,kSi) pairs
		));
	print(tally(L)); --tally the types
	I
	)

--time test(ZZ/5,2,5,3)


--For monomial ideals generated at same degree
--checks their kSS
test = method()
test(Ring,ZZ,ZZ,ZZ) := (k,n,d,extragens) -> (
	R := ZZ;
	L := {}; --initialize empty list
	K := null;
	apply(100, i -> (
		R = randomArtinDegreeMonomialRing(k,n,d,extragens); --random ring
		K = kSS(R,10); --compute kSi up to some bound
		L = L | {K}; --make a list of the kSS lists
		));
	tally(L) --tally the types
	)

--time test(ZZ/5,3,2,1)


-*
Takes a ring S, num of gens n, a degree d, and num of
extra generators extragens
*-
test = method()
test(Ring, ZZ, ZZ, ZZ) := (S,n,d,extragens) -> (
	R         := S[vars(0..(n-1))]; --Polynomial Ring
	monos     := monomialList(R,d);
	standard  := {};
	dataList  := {};
	I         := null; --blank ideal
	K         := null; --blank list
	counter   := 0; --counter for debugging purposes
	for i from 0 to (n-1) do (
		standard = standard | {R_i^d};
		);
	combo := toList(0..(extragens - 1));
	while combo#(-1) < #monos do (
		counter = counter + 1; --increment
		--print counter; --display count
		if isReducedList(monos_combo) then (
			print counter;
			I = ideal(standard | monos_combo); --current ideal
			K = kSi(R/I,10); --compute kSS
			dataList = dataList | {K}; --change the second arg in kSS to get loger res
			print I;
			print K;
			);
		combo = nextList(combo);
		);
	tally dataList
	)

time test(ZZ/5,4,3,1)
--try time test(ZZ/5,4,3,2) with the resolution up to length 10 on campus.
--Looks like it will take a while


restart
loadPackage( "SocleSummands", Reload=>true)
    needsPackage "Points";
time test1(ZZ/101,3,3,1)
time test(ZZ/101,3,3,1)
betti F
Rbar**F.dd_3

r=3;n= 6;

F = resOfRandomPoints(3,27);


Rbar = ringOfRandomPoints(2,6)
kSS(Rbar,8)

netList apply(toList(4..9), p->(p,resOfRandomPoints (binomial(p,2)+p//2)))

netList select(apply(toList(4..50), n -> (
	Rbar = ringOfRandomPoints(3,n);
        mat = (res ideal Rbar).dd_3;
    (degree Rbar,betti mat, kSS(Rbar,3)))), t -> t_2 ==t_2)
toList(4..50)

select({1,2,3}, i-> i==3)

n = 12
time I = ideal (Rbar = ringOfRandomPoints(3,n))
betti res I
time kSS(Rbar,7) -- {4,6} ; last matrix of res I has all quadratic terms. Note: missing r+2 is new.

betti res oo
codim I

Rbar
random(1,Rbar)

----
restart
load  "SocleSummands.m2"
debug SocleSummands
S = ZZ/101[a,b,MonomialOrder => {Position => Down}]

m = random(S^{0,1}, S^{-1,-2,-2})
m' = random(S^{0,1}, S^{-2,-2,-2})
I = minors(2,m)
I' = minors(2, m')
R = S/I
R' = S/I'
kSS(R, 8)
kSS2(R,8)
kSS(R', 8)
kSS2(R', 8)
use R'
F = res(coker vars R', LengthLimit=>8)
F.dd_3

--------------
restart
loadPackage "SocleSummands"
debug SocleSummands
installPackage "Points"
S = ZZ/32003[a,b,c,d]
R = S/ideal"a3,b3,c3,d3,cd2"

S = ZZ/32003[a,b,c]
setRandomSeed 0
elapsedTime scan(1, i-> (
	R = S/(ideal apply (gens S, x->x^3) + ideal random(S^1,S^{-3}));
	<<kSS(R,8)<<flush<<endl;
	))
setRandomSeed 0
elapsedTime scan(1, i-> (
	R = S/(ideal apply (gens S, x->x^3) + ideal random(S^1,S^{-3}));
	<<kSS2(R,8)<<flush<<endl;
	))
setRandomSeed 0
elapsedTime scan(1, i-> (
	R = S/(ideal apply (gens S, x->x^3) + ideal random(S^1,S^{-3}));
	<<kSS3(R,8)<<flush<<endl;
	))

time kSS (R,8)
time kSS2 (R,8)
time kSS (R',8)
time kSS2 (R',8)

--Note: these give different answers! something's wrong.
S = ZZ/32003[a,b,c]
R = S/(ideal apply (gens S, x->x^3) + ideal"a2b")
kSS(R,8)
kSS2(R,8)
kSS3(R,8)

use S
R = S/(ideal apply (gens S, x->x^3) + ideal"abc")
kSS(R,8)
kSS2(R,8)
kSS3(R,8)

F = res coker vars R
m = F.dd_3
socleSummand image F.dd_2
kSS (R,10)
kSi (R,8)
socleSummand2 res (coker vars R, LengthLimit=>8)
socleSummand3 res (coker vars R, LengthLimit=>8)
F = res (coker vars R, LengthLimit=>8)
F.dd_2
vars R
trim ideal generalLinearRow F.dd_2
numgens oo

-------
restart
loadPackage "SocleSummands"
debug SocleSummands
S = ZZ/5[a,b,c]
I = ideal"a3,b3,c3"
R = S/I^3
maxres = 7
elapsedTime kSS2(R, maxres)
elapsedTime kSS(R, 8)

----------

///
installPackage "Points"

--The following code leads to the conjecture that the conductor of the 
--homogeneous coordinate ring of n general points in P^r is always a power
--of the maximal ideal. 

--if binomial(r+k-1,r)<n<= binomial(r+k,r) then the conductor is mm^k

--If this is true, then every indecomposable Ulrich module
--over this ring is the factor ring corresponding to a single point.

--Note: for points on the RNC in P^3, the conductor is not a power of mm starting
--with 8 points.

restart
needsPackage "Points"

S = ZZ/32003[a,b,c,d]
R = S/monomialCurveIdeal(S,{1,2,3})
hilbertFunction(2,R)
mm = ideal vars R
I = trim(mm^4+ideal a^5)
isPowerOfMaxIdeal I

pointsOnRNC = method()
pointsOnRNC (Ring,ZZ) := Matrix => (S,n) ->(
    r := numgens S - 1;
    J := apply(n, j-> random(1000));
    M = map(S^(r+1), S^n, (i,j) -> (J_j)^i)
    )

    

isPowerOfMaxIdeal = method()
isPowerOfMaxIdeal Ideal := Boolean => I -> (
    R := ring I;
    degs := I_*/degree;
    e := degs_0;
    (all(degs, d -> d == e) and 
	numgens I == numcols basis(e,ring I)
    )
)

conductorOfPoints = method()
conductorOfPoints Matrix := Ideal => PM -> (
    Rpoints := ring PM/points PM;
    n := numcols PM;
    L := toList(0..n-1);
    partialLists := apply(L, i -> drop(L,-(n-i))|drop(L,i+1));
    C := intersect apply(n, i-> points(PM_{i})+points PM_(partialLists_i));
    trim sub(C,Rpoints)
    )

S = ZZ/32003[a,b,c,d]
R = S/monomialCurveIdeal(S,{1,2,3})
n = 10

scan(30,n->(
	PM = pointsOnRNC(R,n+2);
	if n+2 == degree points PM then
	<<isPowerOfMaxIdeal conductorOfPoints PM)<<endl)

P = pointsOnRNC(R,20)
ring P
degree points P
ring conductorOfPoints P
isPowerOfMaxIdeal conductorOfPoints P

points randomPointsMat(R,4)
degree oo


loadPackage("SocleSummands", Reload=>true)
needsPackage "Points"
S = kk[a,b,c]

scan(30,i-> (R = S/(points pointsOnRNC(S,i+4));
<<degree R<<"  "<<kSS (R,5)<<endl;
))

scan(15,i-> (R = S/(points randomPointsMat(S,i+4));
<<betti res ideal R<<endl;
<<degree R<<"  "<<kSS (R,5)<<endl;)
)

S = kk[a,b,c,d]
scan(21,i-> (R = S/(points pointsOnRNC(S,i+4));
	<<time betti res coker vars R<<endl;
	))


---
restart
needsPackage "SocleSummands"
needsPackage "DGAlgebras"

S = kk[a,b,c]
n = 2
i0 = ideal apply(gens S, x->x^n)
i = i0+ideal"ab"
i = i
m = 2
im = ideal apply(i_*, g->g^m)

kSi (S/im, 8)
isGolod (S/im)

-------------------------
--Golod cases in 3 variables
restart
needsPackage "DGAlgebras"
needsPackage "SocleSummands"
--Long's test for Golod in monomial ideals in 3 vars.
--(1) [I : x1] · [I : (x2, x3)] ⊆ I for all permutations {x1, x2, x3} of {x, y, z}.
--(2) [I : x1] · [I : x2] ⊆ x3[I : (x1, x2)] + I for all permutations {x1, x2, x3} of {x, y, z}.

///
S = ZZ/5[a,b,c]
I = (ideal random(S^1, S^{-2,-2,-1}))*(ideal random(S^1, S^{-2,-1,-1}))
--I = (ideal random(S^1, S^{-2,-2,-1}))
time isGolod (S/I)
time golodTests I

use S
d = 3
F = random(d,S)
G = random(d,S)
I     = ideal(a^2, b^2,a*F+b*G)

G =  ideal(a^2,b^2, a*c,b*c, c^2-a*b)
G =  ideal(a^2,b^2,c^2)
J = apply(II1, I1 -> I1*G);

II1 = apply(10, i-> ideal random(S^1,S^{-2,-3,-4,-5}));
II2 = apply(10, i-> ideal random(S^1,S^{2:-2,2:-4}));
J = flatten apply(II1, I1 -> apply(II2, I2 -> I1*I2));

apply(#J, n -> kSS (S/J_n,5))
apply(#J, n -> minimalBetti J_n)
kSS(S/J_0,5)

scan(#J, n -> <<n<<"  "<<isGolod(S/J_n)<<"  "<<golodTests J_n<<endl)
scan(#II1, n -> <<n<<"  "<<isGolod (S/J_n)<<"  "<<golodTests J_n<<endl)

golodTests(ideal(a^2,b^2, a*c,b*c, c^2-a*b))
use S

kSS(S/G^2,10)
betti res G
R = S/G^2


------------------------------
needsPackage "SocleSummands"
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

S = ZZ/101[a..c]
L = orbitRepresentatives(S,{3,3,3,3,3,3})
L1 = L/(I-> S/(trim(I + ideal"a3,b3,c3")))
time unique(L1/(R->kSS(R,8)))

L = orbitRepresentatives(S,{3,4,4,5});
L1 = unique (L, I-> trim(I + ideal"a3,b3,c3"));
#L
#L1_0

L2 = L1/(I-> S/I)
#L2
time unique(L2/(R->kSS(R,7)))

S = ZZ/101[a..d]
L = orbitRepresentatives(S,{3,3,4});
L1 = unique apply(L, I-> trim(I + ideal"a3,b3,c3,d3"));
#L1
L2 = L1/(I-> S/I);
time tally(L2/(R->kSS(R,7)))

time unique(L2/(R->kSS(R,7)))
time unique(L2/(R->kSS(R,8)))


--generic example
restart
loadPackage "SocleSummands"
n = 5
S = ZZ/101[x_1..x_n]
d = 3
q = sum(n-1,i->S_i*S_(i+1)) -- works with n = 4, d=2; n= 5, d+3--that is, gives kSS == {n-1,n,...}
--q' = sum(n, i->S_i*S_((i+1)%n)) doesn't work with n= 4!
q = S_0*S_1+S_2*S_3 -- works with n=4, d=2
q = S_0*S_1+S_2*S_3 + S_4*S_5 --- n= 6,d=2 gives {3,4,5...}
I = ideal q + ideal apply(n,i->S_i^d)
kSS(S/I,n+2)

kSS(S/I,2*n-2) == {n-1,n,n+1..2*n-3}

K = koszul vars S
d = 2
e = gens K_2*random(K_2, S^{-2-d})
I1 = ideal (K.dd_2*e)
I2 = ideal(e*K.dd_1)
I = trim(I1+I2)
R = S/I
dim R

primaryDecomposition I
kSS(R,3)
sub(e,R)
syz(vars R)
vars R*e

restart
needsPackage "DGAlgebras"
n = 3
kk = ZZ/101
S = kk[x_1..x_n]
d = 3

I = ideal apply(gens S,y->y^d)
I = ideal(random(S^1,S^{n:-d}))
betti res I
m = ideal vars S
e = (n//2)*d+1
J = truncate(e,I)
J' = truncate(e-1,I)
--J' is not Golod (nontrivial homology product); but, conjecture:J is Golod.
betti res J
R = S/J
A = acyclicClosure(R, EndDegree => 0)
isGolod R
isHomologyAlgebraTrivial A

betti res (m*I)
betti res I

restart
needsPackage "DGAlgebras"
kk = ZZ/101
S = kk[a..e]
I = ideal"ab2, cd2, e3,abcd, d2e2, b2e2, ace, b2d2e"
m = ideal vars S
R = S/(m*I)
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/(m*I))
isHomologyAlgebraTrivial A

betti res I
betti res truncate(4,I)
betti res (m*I)

--random examples quadrics and m^3
restart
loadPackage("SocleSummands", Reload => true)
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

isBurch = R -> (
    if not basis(3,R) == 0 then error ("only works for m^3 = 0 rings");
    A := res coker presentation R;
    ty := rank A_4;
    F = res (coker vars R, LengthLimit => 2);
    rank F_2 > (numgens R)^2 - ty
    )

isBurch1 = R -> (
    I = ideal presentation R;
    m = ideal vars ring I;
    0 != (gens (I:m)) % (m*I:m)
	)
    

S = ZZ/3[a..d]
mm = (ideal vars S)
cubes = ideal apply(numgens S, i->  S_i^3)

flatten for m from 5 to 9 list(
for i from 1 to 10 list(
I = trim(ideal random(S^1, S^{m:-2})+cubes);
R = S/I;
ans:= 0;
print isBurch1 R;
if not isGolod R and not isBurch R then(
ans = kSS(R,8);
<<ans<<endl;
ans)))

isGolod R
isBurch R
kSS(R,8)
kSS(R,8)
isBurch R
----------------------
---
restart
needsPackage "DGAlgebras"
n = 3
kk = ZZ/101
S = kk[x_1..x_n]
d = 3
L = {2,4,7}
I = ideal apply(#L,i-> S_i^(L_i))
m = ideal vars S
betti res I
betti res (m*I)

R = S/(m*I)
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/(m*I))
isHomologyAlgebraTrivial A


--Roos example:
restart
needsPackage "DGAlgebras"
kk = ZZ/101
S = kk[a..d]
I = ideal(a^3, a^2*b, (a+d)*(c^2+b^2), c*d^2, d^3)
m = ideal vars S
R = S/(I)
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/(I))
isGolod (S/(m*I))
isHomologyAlgebraTrivial A


---
--Roos example: Claimed to be non-Golod with trivial homology algebra.
restart
needsPackage "DGAlgebras"
needsPackage "SocleSummands"
kk = ZZ/101
--kk = QQ
S = kk[x,y,z,u]
I = ideal(u^3, x*y^2, (x+y)*z^2, x^2*u+z*u^2, y*y*u+x*z*u, y^2*z+y*z^2) -- has the betti nums as in Roos
--I = ideal(u^3, x*y^2, (x+y)*z^2, x^2*u+z*u^2, z^3+x*z*u, y^2*z+y*z^2) -- has the betti nums in

betti (A = res I)
A.dd_1
A.dd_2
R = S/I
A = acyclicClosure(R, EndDegree => 0)
isGolod (S/I)
isHomologyAlgebraTrivial A

betti res( coker vars R, LengthLimit =>7)

T = QQ[t]
((1+t)^4)*sum(10, i-> (6*t^2+12*t^3+9*t^4+2*t^5)^i)

kSS(R,5)
---Roos "exceptional" examples from 
--Homological properties of the homology algebra of the Koszul complex of a local ring:Examples and questions
restart
needsPackage "SocleSummands"
kk = ZZ/101

S = kk[x,y,z,u]
I = ideal"x2, xy, xz + u2, xu, y2 + z2, zu"
betti res I
R = S/I
kSS(R,7) =={0}

use S
I = ideal "xz + u2, xy, xu, x2, y2 +z2,zu,yz"
betti res I
R = S/I
kSS(R, 7) == {0}

use S
I = ideal"x2 + yz + u2,xy,zu,z2,xz + yu,xu"
betti res I
R = S/I
--kSS(R,9) == {5,6,7,8} --this is fairly slow.

-----Roos examples!
needsPackage"DGAlgebras"
kk = ZZ/101
S = kk[u,x,y,z]
use S
II = {I1 = ideal(u^3, x*y^2, (x+y)*z^2, x^2*u+z*u^2, y^2*u+x*z*u, y^2*z+y*z^2),
-- kSS(R,7)={3,4,5}
I2= ideal"u3,y2u2,z2u2,y2zu,xy2,(x+y)z2,x2yz,x2u",
-- kSS(R,7)={0}
I3=ideal"u3,y2u2,z2u2,y2zu,xy2,(x+y)z2,x2yz,x2u+xyu",
-- kSS(R,7)={0}
I4=ideal"u3, y2u2, z2u2, y2zu, xy2, (x + y)z2, x2yz, x2u + zu2",
-- kSS(R,7)={0}
-- the next two has only cubics
I5=ideal"u3, xy2,(x+y)z2, x2u + zu2, y2u + xzu, y2z + yz2",

-- kSS(R,7)={3,4,5}
I6= ideal"u3, x2u, xz2 + yz2, xy2, x2y+y2u, y2z+z2u"
-- kSS(R,7)={4,5,6}
}
netList for I in II list(
H = acyclicClosure(S/I, EndDegree => 0);
isGolod(S/I), isHomologyAlgebraTrivial(H))
elapsedTime kSS(S/I6, 8)

restart
needsPackage "SocleSummands"
needsPackage "DGAlgebras"
kk = ZZ/101
S = kk[x,y,z]
I0 = ideal apply(gens S, x ->x^3)
I = I0+ideal"xyz"
R = S/I
describe S
isGolod(S/I)
describe S
R === S/I
T = torAlgebra(R)
degrees T
describe T
betti res(coker presentation R)
help torAlgebra
betti res I

torAlgebra(S/I0)
--------------------
restart
--needsPackage "StableResidual"
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"
needsPackage "LocalRings"
1,5,5,1
1,5,6,2
1,5,7,3
1,5,8,4
(almost ci, 1,5,5,1)) 
--help MonomialOrbits

generalLink = J -> (
J0 := gens(J)*random(S^(rank source gens J), S^3);
J' := localTrim (ideal(J0):(J))
)

sumBettis = J ->(
    G := localResolution generalLink J;
    sum(4, i-> rank G_i)
    )

dropBettis = J ->(
F := res J;
G := localResolution generalLink J;
sumBettis F> sumBettis G
)

localTrim = J -> ideal (gens J %(M*J))

--Long's test for Golod in monomial ideals in 3 vars; conjecture for all ideals 
--(1) [I : x1] · [I : (x2, x3)] ⊆ I for all permutations {x1, x2, x3} of {x, y, z}.
--(2) [I : x1] · [I : x2] ⊆ x3[I : (x1, x2)] + I for all permutations {x1, x2, x3} of {x, y, z}.

--this shows: no 5 generator Golod monomial ideals (any fin length golod in 3 vars that contains
--x^a,y^b also contains x^(a-1)y^(b-1), and permutations -- two more gens aren't enough.
    

testGolod = I -> (
    P := {{0,1,2},{1,2,0},{2,0,1}};
    G := gens S;    
    t1 := all(3, i-> (
	    g = G_(P_i);
    gens((I:g_0)*(I:ideal(g_1,g_2)))%I == 0
    and
    gens((I:g_0)*(I:g_1))%(g_2*(I:ideal(g_0,g_1))+I) == 0
)))
testGolod ((ideal gens S)^2)    
isGolod(S/I)



S = ZZ/101[x,y,z]
M = ideal vars S
setMaxIdeal ideal vars S

golodExamples = (S, D, E) -> (
I0 := monomialIdeal apply(#gens S, i -> S_i^(D_i));
L := orbitRepresentatives(S, I0, E);
print(#L);
select(L, I-> isGolod (S/I))
)

use S
test = D ->(
I0 = monomialIdeal ideal"x10,y10,z10";
L= orbitRepresentatives(S,I0,D);
--print(#L);
elapsedTime all(L, I -> not testGolod ideal I)) --3 sec

for j from 5 to 10 do <<(j, all(10, i-> test{j,i+2}))<<flush<<endl;


elapsedTime L/(I -> isGolod(S/I)) --23 sec
random 10
golodExamples(S,{10,10,10},{random 10, random 10})
golodExamples(S,{3,3,3},{3,3})

-------
--random examples quadrics and m^3
--installPackage "Posets"
restart
needsPackage "Posets"
loadPackage("SocleSummands", Reload => true)
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

--note: lin7 is defined in SocleSummands.m2 as the list of all (up to permutation) 
--115 linearly presented monomial ideals of degree 7

S = ring lin7_0

addedMonomials = I -> compress((gens I)%ideal (I0 (degree (first (I_*)))_0))
isLinearlyPresented = I -> (
    degs:= last degrees presentation module I;
    1 + (min(I_*/degree))_0 == (max degs)_0
    )
I0 = d-> monomialIdeal trim (monomialIdeal apply(gens S, a -> a^d) + (ideal(x,y))^d+(ideal(x,z))^d+(ideal(y,z))^d)
ex = d -> numcols basis(d,S) - numgens I0 d;
complement Ideal := I -> (
    firstmon = I_0;
    d := (degree firstmon)_0;
    J := basis(d,S);
    monomialIdeal (J % I)
    )
complements = L -> (
    firstmon = L_0_0;
    d := (degree firstmon)_0;
    J := basis (d,S);
    for I in L list monomialIdeal (J % I)
    )
lin7_0
I' =  complement (lin7_0)
sort flatten ((I'_*)/exponents)
lin7Comp = complements lin7;
lin7Comp/numgens

d = 7
ex 7
J = ideal basis (d,S)
elapsedTime L = for i from 0 to ex(d) list orbitRepresentatives(S,I0(d),toList(i:d)); --d=7: 898 sec.

netList (lin = select(flatten for i from 0 to ex list for j from 0 to #L_i-1 list (
    I := L_i_j;
    (isLinearlyPresented I, addedMonomials I)
    ), p -> p_0===true))

lcmLattice monomialIdeal (lin_0_1)
texPoset oo
netList for i from 0 to ex(d) list for j from 0 to #L_i-1 list minimalBetti ideal L_i_j

L7 = L;

elapsedTime lin7 = select(flatten L7, I -> isLinearlyPresented I); --5 SEC
#lin7
elapsedTime lin7Comps = complements lin7;
I = complement lin7_0
value I

restart
needsPackage "Posets"
--needsPackage "Graphics"
--viewHelp Graphics

I = complement lin7_10

point  = ell -> (
    --ell is a list of 3 integers;
    --output is the coordinate point
    d := 1.0*sum ell;
    (d-1.0*ell_0, 1.0*ell_1+.5*ell_2)
    )
    
plotPoints = method()

outerTriangle = d -> (
    --first coord is horizontal offset from lower edge.
    --second coord is vertical offset from center.
    L := flatten apply(d+1, j-> apply(j+1, 
	         i->( -.5*j + 1.0*i,d-j)));
    <<"\\begin{center}"<<endl;
    <<"\\begin{tikzpicture}[scale=1, vertices/.style={draw, fill=red, circle, inner sep=2pt}]"<<endl;
    for i from 0 to #L-1 do(
    <<"\\node [vertices] ("<<i<<") at "<<L_i<<"{$L_i$};"<<endl;
    );
    <<"\\end{tikzpicture}"<<endl;
    <<"\\end{center}"<<endl;
    )
d = 4
outerTriangle d

exponents Ideal := I -> flatten apply(I_*, m -> exponents m)

plotPoints Ideal := I ->plotPoints exponents I

plotPoints List := L->(
    --L is a list of triples {i,j,k} such as that produced by plotPoints Ideal.
<<"\\begin{center}"<<endl;
<<"\\begin{tikzpicture}[scale=1, vertices/.style={draw, fill=black, circle, inner sep=2pt}]"<<endl;
for i from 0 to #L-1 do(
    <<"\\node [vertices] ("<<i<<") at "<<point (L_i)<<"{};"<<endl;
    );
<<"\\end{tikzpicture}"<<endl;
<<"\\end{center}"<<endl;
    )
outerTriangle 7


exponents I

--4-var experiments
restart
needsPackage "Posets"
loadPackage("SocleSummands", Reload => true)
needsPackage "MonomialOrbits"
needsPackage "DGAlgebras"

kk= ZZ/101
S = kk[x,y,u,v]
d = 8
B = basis(d,S);
expo = {{2,2,2,2},{3,1,2,2},{2,2,3,1}}--lin pres
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I

expo = {{2,2,2,2},{3,1,2,2}}
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I -- lin pres!
--so the fact that the monomials omitted lie in a plane (exp m)_3 ==0
--does not make it "like" the 3 -var case. Not so surprising
--that a "blocker in 3 vars" might not block in 4 vars.

kk= ZZ/101
S = kk[x,y,u]
d = 8
B = basis(d,S);
(e1,e2) = ({2,2,2},{3,1,2})
dist({2,2,2},{3,1,2})
expo = {{2,2,2},{3,1,2}} -- not lin pres
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I

restart
loadPackage("SocleSummands", Reload => true)
kk= ZZ/101
S = kk[x,y,u,v]
d = 8
B = basis(d,S);
expo = {{2,2,2,2},{3,1,2,2},{2,2,3,1}}--lin pres
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I
F = dual res I
betti F

expo = {{2,2,2,2},{3,1,2,2}}
mons = ideal flatten apply(expo, ee -> product(#ee, i->S_i^(ee_i)))
I = ideal(B%mons)
minimalBetti I -- lin pres!

F = res I
G = res ideal B
e = extend( G,F,map(G_0,F_0,1));
prune coker dual e_4

---
S = ZZ/101[a,b,c,d,e]
I = ideal"ac,ad,bd,be,ce"
betti res I
F = res I
F.dd
(F.dd_2)_{4,3,2,1,0}

restart
loadPackage "SocleSummands"
S = ZZ/101[a,b,c]

for i from 1 to 100 do (
    I := ideal random (S^1, S^{7:-3});
    if isLinearlyPresented I then <<endl<<flush<<minimalBetti ideal random (S^1, S^{7:-3})<<endl;
    )

-------Example of a non-burch ks_3
restart
load "SocleSummands.m2"
S = ZZ/101[x,y]
I =ideal"x4,x2y2,y4"
I =ideal"x8,x4y4,y8"
I = ideal"x3, x2y2, y5"
R=S/I
kSS(R,8)
M = coker random(R^2,R^{4:-2})
MSi(M,8)
F = res coker vars R
hasSocleSummand1 image F.dd_3
b = 15
F = res(coker vars R, LengthLimit=>b)
L = for i from 1 to b list (ns=numcols socleSummands coker F.dd_i;
    <<ns<<endl;
    ns)
for i from 1 to b-1 list L_i - L_(i-1)

kSS(R,6)
k=coker vars R
F = res(k, LengthLimit =>5)
F.dd_4
C = image (F.dd_3)_{0..5}
MSi(C,7)
F.dd_3
F = res(C, LengthLimit =>8)
socleSummands image F.dd_2
socleSummands coker F.dd_3
D = image((F.dd_2)_{0..17})
MSi(D,7)
hasSocleSummand D
for i from 1 to 8 list hasSocleSummand coker F.dd_i
M = image F.dd_2
prune socle M

N = coker random(R^2, R^{4:-2})
elapsedTime MSi(N, 10)
F = res(N, LengthLimit => 8)
elapsedTime for i from 1 to 11 list hasSocleSummand image F.dd_i
use S
mmS = ideal gens S
I = ideal"x2,y2"
J =trim( mmS*I + ideal"x2, xy")
J = ideal"x2,xy,y4"
(J:mmS)/J
presentation(I/J)
B = S/J
A = S/I
AB = map(A,B)
presentation (module ideal vars A)
M = image oo
socleSummands module ideal gens B
M' = image (vars B)_{1}
MSi(M', 6)
F = res M'
F.dd_1
M'' = image (vars A)_{1}
MSi(M'', 6)
F = res(coker vars B, LengthLimit=>10)
L = for i from 1 to 10 list numcols socleSummands coker F.dd_i
for i from 1 to 9 list L_i - L_(i-1)
 ----------
 
--What is kSS for n random points in P^2 when they are not burch?
restart
loadPackage "SocleSummands"
loadPackage "Points"
p = 7
n = binomial(p,2) + p//2
ringOfRandomPoints(11,2)

test = n->(
IS = points randomPointsMat(S,n);
I = TS IS;
<< n << "  " <<kSS(T/I, 7)<<endl)
for n from 4 to 40 do test n





test n
n=5
kk = ZZ/32003
--kk = GF(32)

restart
load "SocleSummands.m2"
kk = ZZ/101
S = kk[a, x,y,z]
I = ideal(a^3*x, y^4,z^4, a^3*y - z^3*y, a^3*z - y^3*x)
R = S/I
mmR = ideal gens R
J = ideal(x+y)+mmR^2

MSi(module J,8)


needsPackage SocleSummands
S = kk[x,y,z]
kSS (S/I, 9)
res I
mm
T = kk[x,y]
I = ideal(x^3, x^2*y^2, y^3) 
I = ideal(x^4, x^2*y^2, y^4)
R = T/I
mmR = ideal gens R
J = ideal(x+y)+mmR^2
MSi(module J, 8)

TS = map(T,S,random(T^1,T^{3:-1}));
matrix TS




needsPackage "Points"
S
n=11
IS = points randomPointsMat(S,n);
R = S/IS
mmR^3 = ideal gens R

weakLin(module (mmR^2), 7)

 n-1 then << toString I <<endl<<endl;
    )

n = 3
S = ZZ/101[x_0..x_(n-1)]
T = ZZ/101[x_0..x_1]
TS = map(T,S,random(T^1, T^{n:-1}))
mm = monomialIdeal gens S
d = 5
ze = monomialIdeal(0_S)
L = orbitRepresentatives(S,ideal ze, mm^10, 6);
elapsedTime L = orbitRepresentatives(S,toList(4:4)|toList(7:5));

#L

I = (mm^[2])^2
I = L_2
J = (I*mm):(I:mm);
numcols compress ((gens J) % mm^2)
elapsedTime (L/test)

I = monomialIdeal(x_0^4,x_0^3*x_1,x_0^2*x_1^2,x_1^5,x_1^3*x_2,x_0^3*x_2^2,x_0*x_1^2*x_2^2,x_0^2*x_2^3,x_0*x_1*x_2^3,x_1*x_2^4,x_2^5)
R = S/I
use R
S

data I    

I = ideal"x2,xy,y2"    
data I

I = ideal"x2, xy2, y3"
data I

S = kk[a,b]
I=ideal(a^4, a^2*b, b^2)
data I


loadPackage "SocleSummands"
Lburch3

I = Lburch1_5
S = ring I
JS = (I*mm) : (I:mm)
R = S/I
J = sub(JS,R)
weakLin(J, 7)

MSi(module J,5), 


installPackage "CorrespondenceScrolls"
viewHelp carpet
I = carpet{2,7}
betti res I
hilbertSeries((ring I)^1/I, Reduce =>true)
S = ring I
vars S
kk = coefficientRing S
T = kk[y_0..y_6]
TS = map(T,S,random(T^1, T^{11:-1}))
J = TS I
hilbertSeries(T^1/J, Reduce =>true)

restart
kk = ZZ/101
m=3
n=5
d = 2
S = kk[x_(1,1)..x_(m,n)]
M = genericMatrix(S,x_(1,1),m,n)
J = minors(d,M)
T = kk[y_1..y_((m-d+1)*(n-d+1))]
TS = map(T,S,random(T^1,T^{m*n:-1}))
I  = TS J
mm = ideal gens T
(mm*I) :(I:mm)

------
restart
loadPackage "SocleSummands"
needsPackage "Points"
kk = ZZ/101
S = kk[x,y,z]
I = points randomPointsMat(S, 5)
M = syz gens I
R = S/I
N = sub(M,R)
J = coker transpose (N_{0})
J = coker random(R^1, R^{2:-2})
J = coker matrix"x2,y2"
weakLin(J, 7)
weakLin(J, 7,2)
trim (ideal gens R)^2



F = res J
sy = F.dd_3
sy = coker F.dd_4
prune(H =  Hom(sy,sy))
ho = apply (numgens H, i ->homomorphism H_{i})
T = kk[a_0..a_3]
hoT = apply(ho, A -> sub(matrix A, T))
HoT = sum(4, i-> hoT_i*a_i)
codim minors(2,HoT)
degree minors(2, HoT)
F.dd

vars S

restart
loadPackage "SocleSummands"
needsPackage "AInfinity"
kk = ZZ/101
U = kk[a,b]
m = ideal gens U
I = m^[3]*a*b
I = (ideal"a2,b2")^2
I = (ideal"a3,b3")^2
R = U/I
res coker vars R
mR = sub(m, R)
weakLin(coker random(R^1, R^{-3}),7)
--!
betti res(coker random(R^1, R^{-3}), LengthLimit => 7)
weakLin(module ideal "a3+b3",7)
MSi(module ideal "a3+b3", 7)
MSi(module (mR^2),7)
weakLin(module (mR^2),5)
betti burkeResolution(module ideal "a3+b3", 7)
betti (F = burkeResolution(coker random(R^1, R^{-3}), 7))
picture F
use U
betti res I
RU = map(R,U)
betti (G = res (pushForward(RU, coker random(R^1, R^{-3}))))
G.dd_2

------
needsPackage "MonomialOrbits"
needsPackage "SocleSummands"
kk = ZZ/101
S = kk[a,b,c]
m = ideal gens S
L = orbitRepresentatives(S, m^[5], {4})
netList apply(L, I -> weakLin (S^1/I,5))
I = last L
weakLin(S^1/(last L), 5)
ids = apply(L,I -> (
	F := res I;
	{ideal F.dd_1, ideal F.dd_2, ideal F.dd_3}
	));
for i from 0 to #ids-1 do(
    (gens iid_i_0+ids_i_1  

R = S/m^[3]
J = ideal(a^2*b^2)
weakLin(module J, 7)
m*L_0 : (L_0:m)

-----------
restart
needsPackage "MonomialOrbits"
needsPackage "SocleSummands"
kk = ZZ/101
S = kk[a..f]
m = ideal gens S
L = orbitRepresentatives(S, m^[5], {4})
I = (m^[3])^2
R = S/I
f = a+b+c+d
weakLin (module ideal (f^2),3)
--!
S = kk[x_(0,0)..x_(2,3)]
A = genericMatrix(S, x_(0,0),3,4)
I = minors(3, A)
R = S/I
N = coker random(R^1,R^{-1,-1,-2})
weakLin(N, 4)

restart
loadPackage ("SocleSummands", Reload =>true)
needsPackage "MonomialOrbits"
kk = ZZ/101
c = 2
d = 2

S = kk[x_0..x_(c-1)]
m = ideal gens S
squareFree = (S,d) ->(
    y := symbol y;
    T = kk[y_0..y_(numgens S -1), SkewCommutative => true];
    (map(S,T,vars S))(ideal basis (d,T))
    )
n = squareFree(S,d)
orbitRepresentatives(S,m^[d],m^d,1)


R = S/m^5
F = res coker vars R
prune socle image oo.dd_2
socleSummands coker F.dd_3


--canonical blowup -- when Gorenstein?
installPackage "ReesAlgebra"
installPackage "Points"

restart
needsPackage "ReesAlgebra"
needsPackage "Points"
needsPackage "FastMinors"
--viewHelp

--viewHelp omegaPoints
--viewHelp Points
canonicalIdeal = method()
canonicalIdeal Ideal := Ideal => I ->(
    S := ring I;
    R := S/I;
    F := res I;
    omega := coker sub(transpose F.dd_(length F), R);
    H := Hom(omega,R^1);
    deg := max degrees source gens H;
    g := (gens H)*random(source gens H, R^-deg);
    trim sub(ideal g,R) ---canonical ideal of a 1-dim ring.
)

kk=ZZ/101
S = kk[x,y,z]

I = points randomPointsMat(S,6)
I = minors(2, random(S^3, S^{-1,-2}))


degree I
betti gens I
R = S/I
om = canonicalIdeal I
codim om
RI = reesIdeal om

betti gens om
T = kk[x,y,z, vars(0..numgens ring RI -1)]
vars ring RI
omegaRI = ker map(ring RI/RI, T, vars ring RI|vars R)
betti (F =res omegaRI)
length F
1+numgens ring RI == length F
codim ideal F.dd_(length F)
E = Ext^2(omegaRI,ring omegaRI)
betti prune E
dim T
codim minors(2, presentation E)
--not Gorenstein?

codim chooseGoodMinors(1000, 2, presentation E)
betti prune E
betti res prune E

---weakLin of modules over 0-dim Gorenstein.
restart
loadPackage "SocleSummands"
kk = ZZ/30223
S = kk[x_1..x_4]
I = trim ideal fromDual(sum(4, i-> S_i^2))
betti (F = res I)

R = S/I
f = ideal random(R^1, R^{-2})
betti (F = res f)
ent1 F
weakLin((ring f)^1/f, 5)
M = coker transpose F.dd_4
betti (G = res M)
ent1 G
--!

restart
loadPackage "SocleSummands"
kk = ZZ/30223
S = kk[x_1..x_3]
I = trim ideal fromDual(sum(numgens S, i-> S_i^3))
betti (F = res I)

R = S/I
RS = map(R,S)
M = coker random(R^1, R^{-2})
betti (F = res M)
dim M
ent1 F
weakLin(M, 5)

needsPackage "AInfinity"
--viewHelp picture
mR = aInfinity R;
m =  aInfinity M;

G = burkeResolution(M, 6)
netList for i from 3 to 6 list picture G.dd_i
netList (Km = sort select(keys (m = aInfinity(M)) , k -> class k_0 === ZZ and sum k <5))
displayBlocks G.dd_5
betti res I

--- 2 vars
restart
loadPackage "SocleSummands"
kk = ZZ/30223
S = kk[x_1..x_2]
I = trim ideal fromDual(sum(numgens S, i-> S_i^3))
betti (F = res I)

R = S/I
RS = map(R,S)
M0 = coker (R_0^2)
betti (F = res (M0, LengthLimit => 6))
MM = for i from 1 to 5 list coker transpose F.dd_(i);

FF = for i from 0 to 4 list res(MM_i, LengthLimit => 5)
netList for i from 0 to 4 list (i, betti FF_i, ent1 FF_i)
i = 2
(i, betti FF_i, ent1 FF_i)
G2 = burkeResolution (MM_2, 5)
displayBlocks G2.dd_3
picture G2.dd_3

N = coker transpose F.dd_2
Gmin = res (N, LengthLimit =>5)
netList for i from 1 to length Gmin list  ent1 G.dd_i
netList for i from 1 to length Gmin list  G.dd_i
needsPackage "AInfinity"
--viewHelp picture
mR = aInfinity R;
m =  aInfinity N;

H = burkeResolution(N, 6)
netList for i from 3 to 6 list picture H.dd_i
netList (Km = sort select(keys (m = aInfinity(N)) , k -> class k_0 === ZZ and sum k < 7))
apply(Km, k->(k,ent1 m#k))

netList for i from 1 to 4 list displayBlocks H.dd_i
betti res I

-----
--1-dimensional minimal mult rings

restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z] -- Degrees => {3,7,8}]
T = kk[t]
I = ker map(T,S, {t^3, t^7, t^8})
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.
mm = ideal vars R
F = res mm
M = coker syz vars R
M' = coker syz matrix"x2,y,z"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2
(x^2*M) : M
(x^2*M') : M'
isIsomorphic(M'1, M'++M', Homogeneous => false) --gives error in quasi-hom case.
isIsomorphic(M1,M++M, Homogeneous => false)
(isIsomorphic(M,M1,Homogeneous =>false))_0 == false

-------
restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z,Degrees => {3,8,10}]
S = kk[x,y,z] --,Degrees => {3,8,10}]
T = kk[t]
I = ker map(T,S, {t^3, t^8, t^10})
betti res (S^1/I)
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.

M = coker syz vars R
M' = coker syz matrix"x2,y,z"
M'' = coker syz matrix"x3,y,z"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2
M''1 = coker (res M'').dd_2

isIsomorphic(M'1, M'++M', Homogeneous => false) --gives error in quasi-hom case.
isIsomorphic(M1,M++M, Homogeneous => false)
isIsomorphic(M''1,M'++M'', Homogeneous => false)
(isIsomorphic(M,M',Homogeneous =>false))_0 == false

-------

restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z]
T = kk[t]
I = ker map(T,S, {t^3, t^10, t^11})
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.

M = coker syz vars R
M' = coker syz matrix"x2,y,z"
M'' = coker syz matrix"x3,y,z"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2
M''1 = coker (res M'').dd_2

isIsomorphic(M1,M++M, Homogeneous => false)
isIsomorphic(M'1, M'++M', Homogeneous => false) --gives error in quasi-hom case.
isIsomorphic(M''1,M''++M'', Homogeneous => false)
(isIsomorphic(M,M',Homogeneous =>false))_0 == false

-------
--arf: m
restart
loadPackage "SocleSummands"
kk = ZZ/32003
S = kk[x,y,z,u]
T = kk[t]
I = ker map(T,S, {t^4, t^7, t^9, t^10})
--not arf
R = S/I -- two ulrich modules; the syz of one is two copies of the same one.

M = coker syz vars R
M' = coker syz matrix"x2,y,z,u"
M1 = coker (res M).dd_2
M'1 = coker (res M').dd_2

isIsomorphic(M1,M++M++M, Homogeneous => false)
isIsomorphic(M'1, M++M'++M', Homogeneous => false)
isIsomorphic(M''1,M''++M'', Homogeneous => false)
(isIsomorphic(M,M',Homogeneous =>false))_0 == false

-------
restart
loadPackage "SocleSummands"
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"

kk = ZZ/32003
S = kk[a,b,c]

i = ideal "a4,b4,c4, abc"
kSS(S/i, 8) 
i = ideal "a4,b4,c4, ab3, b2c2"
kSS(S/i, 7) 
i = ideal "a4,b4,c4, a2b2, b2c2"
elapsedTime kSS(S/i, 10) 
i = ideal "a4,b4,c4, ab3,bc3"
isGolod(S/i)
use S
elapsedTime kSS(S/i, 8) 
i = ideal "a4,b4,c4, ab3,b3c"
elapsedTime kSS(S/i, 8) 
i = ideal "a4, a2b2, a2c2"
elapsedTime kSS(S/i, 8) 
isGolod (S/i)

use S
i = ideal"a2,b2"*ideal "a2,b2,c2"
elapsedTime kSS(S/i, 8) 
isGolod (S/i)
use S
golodTests i

use S
I = monomialIdeal "a6,b6,c6"
elapsedTime L = orbitRepresentatives(S, I, {5,5,5,5});#L
LG = select(L, J -> golodTests J and 
                    burchIndex J == 0 );#LG
elapsedTime LG456 = select(LG, J -> kSS(S/J, 6) == {4,5});#LG456		
elapsedTime LG345 = select(LG, J -> kSS(S/J, 4) == {3});#LG345

scan(LG345, I ->(
	R = S/I;
	<<hasSocleSummand ker koszul(2, vars R)<<endl
	)
)
I = LG345_0

LGlin = select(LG, J -> linearFormsInRes J =={true,false})
		    ;#LG
for I in LG do <<I<<" "<< kSS(S/I, 6)<< linearFormsInRes I<<endl
elapsedTime LGlist = for I in LG list{I,kSS(S/I,6), linearFormsInRes I}
select(LGlist, ell -> ell_1!={3,4,5})
LGlist_0
kSS(S/LG_55,8)						    -- 

I = ideal"a3,b3,c3"*ideal"a2,b2,c2"
trim I
burchIndex I
isGolod(S/I)
linearFormsInRes I
kSS(S/I,7)

use S
I = ideal"a3,b3,c3,abc"
trim I
burchIndex I
isGolod(S/I)
linearFormsInRes I
kSS(S/I,7)

kSS(S/LG_0, 8)
F = res I
F.dd_3
golodTests I
F = res LG_0
F.dd_2
R = S/LG_0
G = res coker vars R
G.dd_2

viewHelp EagonResolution
use S
I = ideal"a6,b6,c6, a3b2, b3c2, c3a2"
linearFormsInRes I
kSS(S/I, 7) == {3,4,5}
needsPackage "EagonResolution"
R= S/I

F = eagonResolution(R,5)
picture(F.dd_3,Display => "DisplayBlocks")
mapComponent(F.dd_4,(3,{}),(0,{3}))
picture(F, Display => DisplayBlocks)


---------------
restart
S = ZZ/32003[a_0..a_5]
M = genericMatrix(S, a_0, 2,3)
I = minors(2, M)
R = S/I
M = coker matrix{apply(6, i-> R_i^2)}
F = res(M, LengthLimit => 6);
betti F
tally((flatten entries F.dd_6)/degree)
M1 = submatrixByDegrees(F.dd_5,(7,7),(8,8));

M1 = submatrixByDegrees(F.dd_4,6,7);
betti (G = res (coker M1))
N = coker G.dd_2
M2 = matrix{apply (3, i-> R_i)}
G2 = res coker M2
N2 = coker G2.dd_5
isIsomorphic(N2, N)
betti res N
betti res N2
betti res coker random(R^{1,0}, R^{6:-2})

S = ZZ/32003[a,b]
R = S/(ideal vars S)^3
M = random(R^1,R^{-2})
betti res(coker M,LengthLimit => 10)

-----------
needsPackage "DGAlgebras"
loadPackage"SocleSummands"
netList (L = for i from 1 to 30 list (
    I = ideal ringOfRandomPoints(3,i);
    if burchIndex I == 0 then I else continue
    ))
for I in L do (<<kSS((ring I)/I, 6)<<endl<<flush)
betti res L_2
isGolod (ring L_2/L_2)
degree(ring L_2/L_2)
netList (L = for i from 1 to 20 list (
    I = ideal ringOfRandomPoints(3,i);
    kSS(ring I/I, 6)))
R = ringOfRandomPoints(3,13)
kSS(R,6)
use R


----- in 3 vars: does last map linear imply 
--the second Koszul cycles have socleSummand?
restart
loadPackage "SocleSummands"
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"

kk = ZZ/32003
S = kk[a,b,c]
I = monomialIdeal "a3,b3,c3"
elapsedTime L = orbitRepresentatives(S, I, {5,5,5});#L
LG = select(L, J -> (not golodTests J) and 
                    burchIndex J > 1) ;#LG
elapsedTime LG456 = select(LG, J -> kSS(S/J, 6) == {4,5});#LG456		
elapsedTime LG345 = select(LG, J -> kSS(S/J, 4) == {3});#LG345

J = LG_0
B = drop(flatten entries basis R, 1)
netList (Blin = for f in B list linearDimInRes ideal f)
netList (Blin = for f in B list linearFormsInSyzygy(ideal f,3))
linearFormsInSyzygy(ideal B_3,3)
F = res module ideal B_3
LL = select(L, I->last linearFormsInRes I == true);#LL
elapsedTime M = apply(L, I -> (
	t = last linearFormsInRes I;
	R = S/I;
        s = hasSocleSummand ker koszul(2, vars R);
        (t==s)));
all(M, a -> a)
mm = ideal vars S
apply(LL, I-> (F = res I;ideal F.dd_3 == mm)
)
#M

L_0

use S
I = ideal"a3,b3,c3,abc"
R = S/I
hasSocleSummand ker koszul(2, vars R)
F = res (S^1/I)
F.dd_3
KS 
p = extend(KS, F, p0)
R**p_2
socleSummands(ker koszul(2,vars R))


needsPackage "Bruns"
use S
M = matrix"a,b,c,0;0,a,b,c"
M = random(S^{0,-2}, S^{-1,3:-3})
M = random(S^{0,-1}, S^{4:-3})
M = random(S^{2:0,1:-1}, S^{-2,4:-3})
codim ann coker  M
betti M

I = bruns transpose M
betti (res I).dd_3
	R = S/I;
        s = hasSocleSummand ker koszul(2, vars R)
koszul(3, vars R)
prune image socleSummands ker koszul(2, vars R)

prune image (socleSummands ker koszul(2, vars R)%koszul(3, vars R))
betti((0*R^1):(ideal vars R))
M
betti res (S^1/I)
betti res prune I
codim I
kSS(S/I, 6)
isGolod(S/I)
prune ((0*(ker koszul(2, vars R)):(ideal vars R)))
prune socleSummands ker koszul(2, vars R)
KS = koszul(vars S)
F = res (S^1/I)
p0 = map(KS_0,F_0,1)
p = extend(KS, F, p0)
prune image (R**p_2)
F.dd_3
betti F
betti KS
----
use S
I= ideal"a3,b3,c3"*ideal"a2,b2,c2"

use S
I = ideal"a2b,ab2,abc"
F = res I
KS = koszul(vars S)
p0
p=extend(KS,F,p0)
F.dd
kSS(S/I,7)
-----------
--Is Burch equivalent to having a linear first syzygy in 3 vars?
restart
loadPackage "SocleSummands"
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"

kk = ZZ/32003
S = kk[a,b,c]
I = monomialIdeal "a5,b6,c7"
elapsedTime L = orbitRepresentatives(S, I, {6,6,6});#L
LB = select(L,J->(burchIndex J) == 1 and 
                  (linearFormsInSyzygy(J,2))_1 == 3 );#LB
LL = select(LB, I->linearFormsInRes I == {true,false});#LL

LB_0
betti res oo

elapsedTime M = apply(L, I -> (
	t = last linearFormsInRes I;
	R = S/I;
        s = hasSocleSummand ker koszul(2, vars R);
        (t==s)));
all(M, a -> a)

I = ideal(a^2, b^2, c^2, a*b+b*c+c*a)
R = S/I
        s = hasSocleSummand ker koszul(2, vars R)
isGolod R


restart
needsPackage "SocleSummands"
needsPackage "DGAlgebras"
S = ZZ/101[a,b,c]
I = ideal"a3,b3,c3"*ideal(a,b,c)
R = S/I
basis(S^1/I)
numcols oo
M = subsets flatten entries(basis(5,S^1/I));
M = drop(M,1);
Ilist = apply(M,m -> sub(monomialIdeal m, R));#Ilist
elapsedTime L =  for m in Ilist list (
    L := MSiList(R^1/m,6);
    s := all(L_{0..4},t -> not t);
    if s then (<<(m, L)<<endl;m)
    );#L

elapsedTime L =  for m in Ilist list (
    L := linearFormsInSyzygy(m,4);
    if L_2 !=3 then 
    <<(m, L)<<endl;m);#L

linearFormsInSyzygy(ideal"abc",4)
isGolod R
(res(S^1/I)).dd
MSiList(R^1/ideal"abc", 6)

N = coker random(R^2,R^{2:-2})
weakLin(N,5)

elapsedTime MSiList(N, 7)
betti (F = res(N, LengthLimit => 5))

I = ideal"a3,b3,c3"*ideal(a,b,c)
R = S/I
linearFormsInSyzygy(ideal"abc", 4)
F = res module ideal "abc"
betti F
isGolod R

--
restart
loadPackage "SocleSummands"
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"
S = ZZ/101[a,b,c,d]
I = ideal"a3,b3,c3,d3"*ideal(a,b,c,d)
R = S/I
isGolod R
linearFormsInSyzygy (ideal"abcd",5)
B = drop(flatten entries basis R,1);#B
subsets(B,2)
S2B = for s in subsets(B,2) list(
          if numgens (I := trim ideal s) == 1 then continue
	  else I)
      #S2B

for b in S2B do(
	Lb := linearFormsInSyzygy( b, 4);
	if Lb_2 !=4 then <<b <<" "<<Lb<<endl)
MSi(module ideal"abcd", 6)
betti (F = res module ideal"abcd")
trim ideal F.dd_3

I = monomialIdeal "a3,b4,c3,d3"
elapsedTime L = orbitRepresentatives(S, I, {3,3,3,3});#L
LB = select(L,J->(burchIndex J) == 1 and 
                  (linearFormsInSyzygy(J,1))_1 == 4 );#LB
LL = select(LB, I->linearFormsInRes I == {true,false});#LL
LB_0
LB_1
R = S/LB_0
basis R
M = coker random(R^1, R^{-3})
J = ideal"abcd"
betti res J
linearFormsInSyzygy(M,6)
MSi(M,5)
MSi(module J, 5)
M = map(R^2, R^5, matrix"a,b,c,d,0;0,a,b,c,d")
MSi(coker M, 6)
linearFormsInSyzygy(coker M,6)

use S
I = monomialIdeal "a3,b3,c4,d4"
elapsedTime L = orbitRepresentatives(S, I, {5,6,6,6});#L
Mlist = apply(#L, i -> (
	I = ideal select(L_i_*, m -> degree m != {5});
	R = S/I;
	M = coker (R**((gens L_i)%I));
	M
	));
elapsedTime Mb = select(Mlist, M -> burchIndex ideal presentation ring M >0)
#Mb
Mb0 = select (Mb, M-> 4>lin1 syz presentation ring M);
#Mb0


elapsedTime Mb3 = select (
    apply(Mb, M-> lin1 res(M, LengthLimit => 4)),
       li -> li_3 != 4);
#Mb3

lin1 res (Mb_1, LengthLimit => 4)


	burchIndex  li -> li_4 !=4)

LB = select(L,J->(burchIndex J) == 1 and 
                  (linearFormsInSyzygy(J,2))_1 == 3 );#LB
-------------------
--3 var pattern of linearity vs Burch index
restart
uninstallPackage "SocleSummands"
restart
installPackage "SocleSummands"

restart
loadPackage ("SocleSummands", Reload =>true)
lin1 = method()
lin1 Matrix := f -> 
      numgens trim ideal ((gens ideal f) % (ideal vars ring f)^2)
lin1 ChainComplex := F -> apply(length F, i -> lin1 (F.dd_(i+1)))
lin1 Module := M -> lin1 res M
lin1 Ideal := I -> lin1 res module I

    needsPackage "DGAlgebras";

data = method()
data Ideal := I ->(
--g := golodTests I;
--g := isGolod ((ring I)/I);
R := S/I;
ell := lin1 res I;
K := koszul vars R;
n := numgens R;
c = for i from 1 to n-1 list hasSocleSummand ker K.dd_i;
soc := kSS (S/I, 6);
--{g,ell,c,soc}
{ell,c,soc}
)
data List := L -> apply(L, J->data J)
needsPackage "MonomialOrbits"

lessData = method()
lessData Ideal := I ->(
--g := golodTests I;
--g := isGolod ((ring I)/I);
R := S/I;
ell := lin1 res I;
K := koszul vars R;
n := numgens R;
c = hasSocleSummand ker K.dd_(n-1);
{ell,c}
)
lessData List := L -> apply(#L, i->(i,lessData L_i))
needsPackage "MonomialOrbits"


kk = ZZ/32003
S = kk[a,b,c]
use S
I =ideal"a3,b3,c3"*ideal"a2,b2,c2"
elapsedTime data I


use S
I = ideal"a4,b4,c4"
I1 = I + ideal"ab3"
I2 = I + ideal"ab3, cb3"
I3 = I+ ideal"ab3,bc3"
I4 = I+ ideal"a2b2,b2c2"
I5 = I+ideal"abc"

L = {I1,I2,I3,I4,I5}
netList data L
use S
I = monomialIdeal "a4,b4,c4"
I = ideal"a5,b5,c5"
elapsedTime L = orbitRepresentatives(S, I, {5,5,6,6});#L
elapsedTime D = lessData L;
netList select(D, s-> s_1_1 == false)

LG = select(L, J -> golodTests J)
#LG
netList data LG
LB = select(L,J->(burchIndex J) == 1 and 
                  (linearFormsInSyzygy(J,2))_1 == 3 );#LB
LL = select(LB, I->linearFormsInRes I == {true,false});#LL

--4 vars golod
S = kk[a,b,c,d]
I = monomialIdeal"a3,b3,c3,d3"
elapsedTime L = orbitRepresentatives(S, I, {3,3,4});#L
--LG = select(L, J -> isGolod((ring J)/J)
netList data L

use S
I = ideal"a3,b3,c3"*ideal vars S
data I
I = monomialIdeal "a5,b6,c7"


elapsedTime L = orbitRepresentatives(S, I, {6,6,6,6});#L
elapsedTime LB = select(L,J->(burchIndex J) == 2);#LB
elapsedTime tally LB/(J->lin1 module J)


     and 
                  (linearFormsInSyzygy(J,2))_1 == 3 );#LB
LL = select(LB, I->linearFormsInRes I == {true,false});#LL

LB_0
betti res oo

elapsedTime M = apply(L, I -> (
	t = last linearFormsInRes I;
	R = S/I;
        s = hasSocleSummand ker koszul(2, vars R);
        (t==s)));
all(M, a -> a)

----
restart
uninstallPackage "SocleSummands"
restart
installPackage "SocleSummands"

restart
loadPackage ("SocleSummands", Reload =>true)
lin1 = method()
lin1 Matrix := f -> 
      numgens trim ideal ((gens ideal f) % (ideal vars ring f)^2)
lin1 ChainComplex := F -> apply(length F, i -> lin1 (F.dd_(i+1)))
lin1 Module := M -> lin1 res M
lin1 Ideal := I -> lin1 res module I

    needsPackage "DGAlgebras";

lindata = method()
lindata (Ideal, ZZ) := (I,r) ->(
--g := golodTests I;
--g := isGolod ((ring I)/I);
R = ring I;
ell := lin1 res( I,LengthLimit => r);
--K := koszul vars R;
ell,c,soc
)
lindata (List, ZZ) := (L,r) -> apply(L, J->lindata (J,r))
needsPackage "MonomialOrbits"
--lin1 for powers of mm?
lin1
S = kk[a,b]
mm = ideal gens S
I = mm^5
R = S/I
J = ideal random(R^2, R^{-2,-3,-4})
J = (ideal gens R)^2
lin1 res (module ideal"abc")
lin1 res (module ideal random(R^1,R^{-2}))
elapsedTime MSi(R^1/J, 6)
viewHelp
kk = ZZ/32003
S = kk[a,b,c,d]
mm = ideal S_*
use S
R= S/mm^6
I = ideal random(R^1, R^{-3,-3})
I = ideal"ab,cd"
lindata(I, 3)
lin1 res (R^1/(ideal"a2,b2"))
lindata(ideal("a2,b3"),3)
installPackage "AInfinity"
golodBetti(R^1/(ideal"a2,b2"),4)

m =  random(R^2, R^{2:-1})
betti res(coker syz m, LengthLimit => 3)
golodBetti(coker syz m,3)

S = kk[a,b,c]
mm = ideal S_*
I = ideal"a3,b3,c3"*mm
R = S/I
M = R^1/(ideal"abc")
F = res(M, LengthLimit => 5)
MM = apply(length F + 1, i-> coker F.dd_(i+1));

netList apply(4, i->(
	N := coker F.dd_(i+1);
betti res(N, LengthLimit => 3),
golodBetti(coker F.dd_(i+1),3), lin1 F)
)
---
--2 vars: is every2nd syz golod?
--canonical+3points
restart
kk = ZZ/101
P2 = kk[a,b,c]
points3 = ideal"ab,ac,bc"
gpoints = gens points3
Cideal = ideal( gpoints * random(source gpoints, P2^{-4}))
C = P2/Cideal
pointsOnC = sub(points3, C)
Kstar = ideal(a+b+c)
Kstarminus3 = (Kstar*pointsOnC)
ID = Kstarminus3

burchOfDivisor = ID -> (
    --ID is the ideal of a divisor of an ideal on a variety. 
    --Take the dual, represented as an ideal, by multiplication
    --with one of its generators, of degree e
    -- and its gens of degree e
    --as the variables...
C = ring ID;
f = ID_*_0;
e = degree f;
J =ideal(f): ID;
L = (gens J)*map(source gens J,,matrix basis(e,J));
<<numcols L<<endl;
T = kk[z_0..z_(numcols L-1)];
IC = trim ker map(C,T,L);

T' = kk[x_0..x_(numcols L-3)];
phi = map(T',T,random(T'^1, T'^{numcols L:-1}));
Ibar = phi IC;
needsPackage "SocleSummands";
(IC,Ibar, burchIndex Ibar)
)


g= 5
ps = points randomPointsMat(P2,3, AllRandom=>true)
restart
curveOfGenus456withdivisor = (g,ps) ->(
    --return
    --C ring of a plane curve of genus g
    --(represented as a quintic with 6-g nodes)
    -- with equation contained in ideal ps of points E
    --K the canonical divisor of the curve (adjoint through the nodes)
    --ID the ideal of the divisor D = K+E on C.
needsPackage "Points";
P2 = ring ps;
sings = if g == 6 then ideal (1_P2) else points randomPointsMat(P2,6-g);
constraint = gens intersect(ps, sings^2);
IC = ideal(constraint*random(source constraint, P2^{-6}));
RC = P2/IC;
singsC = gens trim sub(sings, RC);
IKC =  ideal(singsC*map(source singsC, , matrix basis(2,ideal singsC))); --d-3
Cps = sub(ps, RC);
(IKC, Cps)
)

--ps = points randomPointsMat(P2,3,AllRandom=>true)
globalSectionsInverse = ID->(
    --ID is the ideal of a divisor
    --finds the global sections of ID^{-1}
    f := ID_*_0;
    e := degree f;
    J = ideal f:ID;
    L = (gens J)*map(source gens J,,matrix basis(e,J))
    )
P2 = ZZ/101[a,b,c]
ps = ideal 1_P2
g = 5
(IKC,Cps) = curveOfGenus456withdivisor(4,ps)
numgens IKC
globalSectionsInverse IKC    
numgens trim ideal oo


---
needsPackage "MonomialOrbits"
needsPackage "SocleSummands"
S = ZZ/101[a..d]
I0 = monomialIdeal apply(gens S, x -> x^3);
J = monomialIdeal ((gens (ideal gens S)^4)%I0);
L = orbitRepresentatives(S,I0,{3,3,3,3,3});#L
--L = orbitRepresentatives(S,I0,J, -6);#L
LnotB= select(L, ell -> burchIndex ell ==0);#LnotB
LB  =  select(L, ell -> burchIndex ell >0);#LB
LG = select(L, ell -> isGolod(S/ell));#LG
LGB = select(LG, ell -> burchIndex ell >0);#LGB
LNB = select(L, ell -> burchIndex ell ==0);#LNB
koS LNB_0
kSi (LNB_30, 6)
isGolod(S/LNB_0)

--A non-Golod, non-burch ideal in S = ZZ/101[a..d]. 4th syz of k has socle summand. 5th?
I = monomialIdeal(a^3,b^3,a^2*b*c,a*b^2*c,a^2*c^2,b^2*c^2,c^3,a^2*b*d,a*b^2*d,a*b*c*d,a*c^2*d,b*c^2*d,a^2*d^2,b^2*d^2,a*c*d^2,b*c*d^2,d^3)
kSi(I)
--needsPackage "DGAlgebras"
I = L_1
kSS (I, 5)
koS I

-----Test the generation of Koszul socle generators ----
restart
needsPackage "SocleSummands"
a  = symbol a
S = ZZ/101[a_1..a_6]
mS = ideal vars S
powers =(S, d) -> ideal apply(gens S, x -> x^d)
powers(S,3)
R = S/(mS^2*powers(S,3));
R = S/((powers(S,3))^2);
socleSummandGeneration R

--elapsedTime socleSummandGeneration R -- 93 sec with 10 vars.


S = ZZ/101[a,b,c]
mS = ideal vars S

I2 = ideal apply( gens S, x ->x^2)
I3 = ideal apply( gens S, x ->x^3)
I6 = ideal apply( gens S, x ->x^6)
I = I2*I3
I = I6+ideal"a3b2,b3c2,c3a2"

R = S/I
kSi (I, 5)
isGolod R

betti res I
n = numgens R
mmR = ideal vars R
K = koszulComplexDGA R
L = toComplex K
A = K.natural
soc = apply(n, i-> gens ((0_R*L_(i+1)):mmR));
netList soc
cyc = apply(n, i-> syz K.dd_(i+1));
netList cyc
image K.dd_3
i = 1 
soc_1%(vars R**cyc_i)
socSummands = apply(n, i-> soc_i%(vars R**cyc_i));
netList socSummands


(dgAlgebraMultMap(K,A_2))_1*socSummands_1

socleSummands ker L.dd_2
prune image socSummands_0

cyc_0*(soc_0//cyc_0)
soc_0%(vars R ** cyc_0)
multSocSummands = (dgAlgebraMultMap(K,A_0))_1*socSummands_0||(dgAlgebraMultMap(K,A_1))_1*socSummands_0
(vars R)
 cyc_1
multSocSummands % (R^{-1}**vars R **cyc_0)
((dgAlgebraMultMap(K,A_1))_1*socSummands_0)%(vars R **cyc_1)

prune image (vars R**syz L.dd_1)
prune image socSummands_0
L_1 == target syz K.dd_1
L_1 ==ambient cyc_0
L_1 == ambient soc_0


--A non-generation example in a non-Golod ring
--with no socle summands in the resolution, 
--though there are in the Koszul cycles.
S = ZZ/101[a,b,c]
mS = ideal gens S
 I1 = ideal (a^4,b^4,c^4, a*b^3)
R = S/I1

R = S/(mS*I1)
koSi R
kSi R
isGolod R
R = S/LB_0
socleSummandGeneration R
kSi(R,4)

isBurch LB_0

 ---
--Investigate persistence of socle summands in Koszul cycles, esp in Golod rings.
--Assume R = S/I Golod.
--If Burch then there should be socle summands from syz_2 on
--if not burch and cyc_3 has a socle summand, does cyc 4?
--we know that cyc_(numgens) is the socle, so we need
--at least 5 variables to see this phenomenon.
restart
needsPackage "SocleSummands"
needsPackage "DGAlgebras"
needsPackage "MonomialOrbits"
persists = ell -> (
    ell' := apply(ell, t-> if t then 1 else 0);
    ell' = drop(ell', 1);
    m := min positions(ell', i-> i== 1);
    sum ell' == length ell'-m
    )

kk= ZZ/101
S = kk[a..e]
I0 = ideal apply(numgens S, i-> S_i^2)
I1 = ideal apply(numgens S, i-> S_i^3)
L = orbitRepresentatives(S, I1, toList(3:2)|{3})
L = orbitRepresentatives(S, I1, toList(3:2))
#L
elapsedTime L1 = (for I in L list kkSi(I^2));
Lt = apply(L1, ell -> persists ell);
positions(Lt, t -> false)


----
loadPackage( "SocleSummands", Reload => true)
kk = ZZ/101
S = kk[x,y,z]
I = ((ideal vars S)*(ideal"x3, y3, z3"))
R = S/I
M = R^1/ideal"xyz"
betti (F = res( M, LengthLimit =>8))
kSS (R, 5)
burchIndex I
MSi(M, 7)
MSi(coker vars R, 4)

---
--look for rings of burch index 2 with a module M whose
--6-th syz has no socle summand (then the 5th and 7th will have these).
restart
loadPackage( "SocleSummands", Reload => true)
kk = ZZ/101
S = kk[a,b]
I = minors(2, matrix"a,b,0;b8,ab7,a8")
socle(S^1/I)
--I = minors(2, A = random(S^2, S^{-1,-4,-3}))
burchIndex I
R = S/I
M = coker random(R^2,R^{-4,-4, -5})
M = coker random(R^2,R^{-4,-4, -5})
M = coker random(R^2,R^{-3,-2})
    MSi(M, 6)
burchIndex ideal ring M
F = res

restart
loadPackage( "SocleSummands", Reload => true)
kk= ZZ/101
S = kk[a..c]
I0 = ideal apply(numgens S, i-> S_i^7)
L = orbitRepresentatives(S, I0,{3,3,3,5})
#L
L2 = select(L, ell -> burchIndex ell >1);#L2
I = L2_4
R = S/I
M = coker random(R^2,R^{-4,-5,-6})
MSi(M,6)

----
restart
loadPackage "SocleSummands"
--loadPackage ("DGAlgebras", Reload =>true)
uninstallPackage "SocleSummands"
restart
installPackage "SocleSummands"
kk = ZZ/101
S = kk[a,b,c]
I = ideal"a6,b6,c6, a3b2, a3c2"
isGolod(S/I)
betti res I
R = S/I
K = koszul(vars R)
use R
kSi (F = res (coker vars R, LengthLimit => 8))
code methods kSi
socleSummands image F.dd_4
socleSummands image F.dd_6
socleSummands image F.dd_8
prune socle (S^1/I)
A = koszulComplexDGA R
keys (A.cache)
AA = A.natural
AA_0*AA_1*AA_2
B = A.cache.homologyAlgebra
C = B#homologyAlgebra
presentation C
mC = ideal (vars C)
trim (mC^3)

hilbertFunction(0,C)
AA = A#UnderlyingAlgebra
A_0
A_0*A_1*A_2
viewHelp DGAlgebras
B = HH A
presentation B
----



restart
needsPackage "SocleSummands"
code methods killCycles
kk = ZZ/101
S = kk[a,b]
R = S/ideal"a2,ab,b2"
K = koszulComplexDGA R
K1 = killCycles(K, EndDegree => 1)
K2 = killCycles(K1, EndDegree => 2)
K3 = killCycles(K2, EndDegree => 3)
K4 = killCycles(K, EndDegree => 4)
K5 = killCycles(K, EndDegree => 5)
toComplex(K5, 5)
toComplex (acyclicClosure(K, EndDegree=>4), 5)

elapsedTime netList (KK = KtoR (R,1,5))
extend(KK_1,KK_0, map(KK_1_0, KK_0_0, 1))
socleSummands(KK_0)
socleSummands(KK_1)
socleSummands (res( coker vars R, LengthLimit => 5), 5)
socleSummands res coker vars R
viewHelp ChainComplex


--define a set of degree -1 maps on the tate resolution; see if they are jointly surjective
--(these are extensions -- not necessarily extensions as derivations; but those do exist,
--and must be homotopic to these.)
restart
S = ZZ/101[a,b]
R = S/ideal"a3,a2b,b3"
F = res(coker vars R, LengthLimit => 8)
G = chainComplex for i from 1 to max (F[1]) list (F[1]).dd_i;
e = apply(numgens G_0, t -> map(F_0, G_0 ,(i,j) -> if j == t then 1_R else 0));
E=apply(numgens G_0, t -> extend(F,G,e_t, Verify => true));
Etot1 = apply(1+length G, t -> apply(numgens G_0, ell -> E_ell#t));
Etot = apply(1+length G, ell  -> 
     map(F_ell, R^(numgens G_0) ** G_ell, fold((a,b) -> a|b, Etot1_ell)));
apply (1+length G, ell -> isSurjective Etot_ell)

restart
needsPackage "SocleSummands"
n=4
S = ZZ/101[x_1..x_n]
product gens S
I = ideal fromDual (product gens S + sum apply(gens S, x ->x^(numgens S)))
betti res I
R= S/I
koszul vars R
socleSummands koszul vars (S/I)
KtoR(R,3, 5)



---CM modules over points
needsPackage "Points"
viewHelp Points
n = 3
S = ZZ/32003[x_0..x_n]
M = randomPointsMat(S,n+1)
I = points M
R = S/I
ipoints = subsets(toList(0..n), n)
ideals = apply(ipoints, s -> ideal((gens R)_s))
MM = ideals/(I ->R^1/I)
omega = Ext^n(S^1/I, S^{-n-1})
(gens R)_(ipoints_0)


S = ZZ/101[x_1..x_4]
I1 = ideal(x_1^2, x_2^2, x_3^2, x_4^2, x_1*x_2, x_2*x_3, x_3*x_4)
I = ideal(x_1^3, x_2^2, x_3^2, x_4^2, x_1*x_2, x_2*x_3, x_3*x_4)
I = I + ideal(x_1^2*x_3, x_1^2*x_4)
numgens I
numgens trim I
trim I
R = S/I
isBurch I
socleSummands res(coker vars R, LengthLimit =>5)
mmS = ideal gens S
gens (mmS^3) % I

I = I1
R = S/I
mmR = ideal vars R
isBurch I -- false
isGolod R -- false
socleSummands res(coker vars R, LengthLimit =>5)
mmS = ideal gens S
gens (mmS^3) % I
M = R^1/(x_1*mmR)
M = coker random(R^1, R^{2:-2})
socleSummands res(M, LengthLimit =>5)


--
