-*
newPackage("MultigradedBGG",
    Version => "1.1",
    Date => "5 June 2023",
    Headline => "Computing with the multigraded BGG correspondence",
    Authors => {
	{Name => "Maya Banks",         	     Email => "mdbanks@wisc.edu",      HomePage => "https://sites.google.com/wisc.edu/mayabanks" }
        {Name => "Michael K. Brown",         Email => "mkb0096@auburn.edu",    HomePage => "http://webhome.auburn.edu/~mkb0096/" }
	{Name => "Daniel Erman",    	     Email => "erman@wisc.edu",        HomePage => "https://people.math.wisc.edu/~erman/" }
	{Name => "Tara Gomes",	    	     Email => "gomes072@umn.edu",      HomePage => "Fill in" }
	{Name => "Prashanth Sridhar",	     Email => "pzs0094@auburn.edu",    HomePage => "https://sites.google.com/view/prashanthsridhar/home" }
	{Name => "Eduardo Torres Davila",    Email => "torre680@umn.edu",      HomePage => "https://etdavila10.github.io/" }
	{Name => "Sasha	Zotine",    	     Email => "18az45@queensu.ca",     HomePage => "https://sites.google.com/view/szotine/home" }
	    },
  DebuggingMode => true
  )


--Not sure about all these exports--sort out later.
exports {
    "addTateData",
    "dualRingToric",
    "RRFunctor",
    "sortedMonomials",
    "positionInKoszulComplex",
    "dualElement",  --keep.
    "addZeroTerms", -- drop.
    "completeToMapOfChainComplexes", -- might use addTerms in an option.  drop "true" option and option altogether?
    "greaterEqual",
    "strictlyGreater",  --drop.
    "degreeTruncation",
    "isChainComplexMap", -- should be in the core of M2
    "degreeTruncationAbove", -- working from the other side, not tested
    "degreeSetup", -- get the degree range from a starting set of degrees
                   -- using c multiplications with variabels in S
    "strictDegreeTruncationAbove",  --drop.
    "concatMatrices",
    "matrixContract",
    "tallyComplex",
    "annHH", --drop
    "dimHH", --drop
    "relevantAnnHH",  --keep
    "beilinsonWindow",
    "factors", -- now superflous, because the last degree component in E
               -- give this information
    "entry", 
    "DMonad", --drop
    "cacheComplexes",
    "cachePhi",
    "diffModToChainComplexMaps",  
    "bigChainMap",
    "RRTable",
    "doubleComplexBM",
    "isIsomorphic",
    "Variable",
    "SkewVariable",
    "unfold"
    }
*-


--path=prepend("../",path)
load "DifferentialModules.m2";
---
---
---
---

--Input:    A pair of matrices (M,N)
--Output:   The effect of contracting M by N. 
matrixContract = method()
matrixContract (Matrix,Matrix) := (M,N) -> (
    S := ring M;
    if M==0 or N==0 then return map(S^(-degrees source N),S^(degrees source M),0);
    assert(rank source M == rank target N); 
    transpose map(S^(-degrees source N), , transpose matrix apply(rank target M,i->apply(rank source N,j->		
           sum(rank source M,k->contract(M_(i,k),N_(k,j) ))
	    )))
    )


dualRingToric = method(
        Options => {
    	    Variable        => getSymbol "x",
	    SkewVariable => getSymbol "e"
	}
    );

dualRingToric PolynomialRing := opts -> S ->(
    --  Input: a polynomial ring OR an exterior algebra
--  Output:  the Koszul dual exterior algebra, but with an additional 
--           ZZ-degree, the ``standard grading'' where elements of \bigwedge^k
--           have degree k.
    kk := coefficientRing S;
    degs := null;
    if isSkewCommutative S == false then(
    	degs = apply(degrees S,d-> (-d)|{-1});
	e := opts.SkewVariable;
    	ee := apply(#gens S, i-> e_i);
    	return kk[ee,Degrees=>degs,SkewCommutative=>true, MonomialOrder => Lex]
	); 
    if isSkewCommutative S == true then(
    	degs = apply(degrees S, d-> drop(-d,-1));
	y := opts.Variable;
    	yy := apply(#gens S, i-> y_i);
    	return kk[yy,Degrees=>degs,SkewCommutative=>false]
	);       
    )

TEST ///
restart
needsPackage "NormalToricVarieties"
load "MultigradedBGG.m2"
S = ring(hirzebruchSurface(2, Variable => y));
E = dualRingToric(S,SkewVariable => f);
SY = dualRingToric(E);
assert(degrees SY == degrees S)
///

toricRR = method();
--Input: (M,LL) M a (multi)-graded S-module.
--             LL a list of degrees.
--Output:  The differenial module RR(M) in degrees from LL,
--         presented as a complex in homological degrees -1, 0 ,1
--         and with the same differential in both spots.
toricRR(Module,List) := (N,LL) ->(
    M = coker presentation N;
    S := ring M;
    if not isCommutative S then error "ring M is not commutative";
    if not S.?exterior then S.exterior = dualRingToric(S);
    E := S.exterior;
    relationsM := presentation M;
    -- this used to say "gens image presentation M"... just in case a bug arises
    f0 := matrix {for d in unique LL list gens image basis(d,M)};
    wEtwist := append(-sum degrees S, -numgens S);
    df0 := apply(degrees source f0, d -> (-d | {0}) + wEtwist);
    SE := S**E;
    --the line below is better for degrees,it overwrites S somehow...
    -- sasha edit: i edited the line below so that it wouldn't overwrite S. dunno if we want this.
    --SE := (coefficientRing S)(monoid[gens S|gens E, Degrees => apply(degrees S,d->d|{0}) | degrees E, SkewCommutative => gens E]);
    tr := sum(dim S, i-> SE_(dim S+i)*SE_i);
    newf0 := sub(f0,SE)*tr;
    relationsMinSE := sub(relationsM,SE);
    newf0 = newf0 % relationsMinSE;
    newg := matrixContract(transpose sub(f0,SE),newf0);
    g' := sub(newg,E);
    if E^df0 == E^0 then chainComplex map(E^0, E^0, 0) else (
    	differentialModule(chainComplex{map(E^df0,E^df0, -g', Degree => degree 1_S | {-1}),map(E^df0,E^df0, -g',  Degree => degree 1_S | {-1})}[1])
    	)
    )

-- if no degree range is inputted, makes a simple choice.
toricRR Module := M -> (
    LL := if length M < infinity then unique flatten degrees basis M
    else join(degrees M,apply(degrees ring M, d -> (degrees M)_0 + d));
    toricRR(M,LL)
    )

TEST ///
restart
loadPackage "NormalToricVarieties"
load "MultigradedBGG.m2"
S = ring hirzebruchSurface 3;
M = coker matrix{{x_0}};
L = {{0,0}, {1,0}, {-3, 1}, {0,1}, {2,0}};
D = toricRR(M,L)
D.dd


LL = {{1,0}, {-3, 1}, {0,1}}
toricRR(N, LL)
LL = {{0,0}, {1,0}}

RM = toricRR(M, LL)
toricRR(M)

S = ring weightedProjectiveSpace {1,1,1,1}
N = coker map(S^1, (S^{-2})^4, matrix{{x_0^2, x_1^2, x_2^2, x_3^2}})
isHomogeneous N
M = coker map(N**(S^{2}) ++ N**(S^{1}), N**(S^{1}) ++ N ++ N ++ N, matrix {{x_0, 0, 0, x_1*x_3}, {0, x_3, x_1, -x_0}})
isHomogeneous M

basis M
LL = {-2,-1, 0,1}
toricRR(M, {-2,-1, 0,1})
oo.dd
X = weightedProjectiveSpace {1,1,2}
S = ring X
M = coker matrix{{x_0, x_1^2, x_2}}
F = res M
F.dd
presentation M
flatten degrees target oo
toricRR(M, {0,1,2,3})
///

TEST ///
restart
load "MultigradedBGG.m2"
kk=ZZ/101
S=kk[x_0,x_1,Degrees=>{1,1}]
D = toricRR(S^1,{0,1})
G = cornerDM({0},D,LengthLimit => 3)
G_0 == minimalPart G
tally degrees G_0
F = resDM(D, LengthLimit => 2)
tally degrees F_0
tally degrees minimalPart F
F.dd_1
assert(D.dd^2 == 0)
assert(isHomogeneous D)
///


TEST ///
restart
load "MultigradedBGG.m2"
kk = ZZ/101
-- ring of hirzebruchSurface 3
S = kk[x_0, x_1, x_2, x_3, Degrees =>{{1,0},{-3,1},{1,0},{0,1}}]
M = S^{{2,-4}}
L = unique apply(((ideal 1_S)_* | (ideal {x_0..x_3})_* | ((ideal {x_0..x_3})^2)_*)/degree, x -> x + {-2,4})
RM = toricRR(M,L)
assert(RM.dd^2 == 0)
assert(isHomogeneous RM)
///

TEST ///
restart
load "MultigradedBGG.m2"
kk = ZZ/101
-- ring of hirzebruchSurface 3
S = kk[x_0, x_1, x_2, x_3, Degrees =>{{1,0},{-3,1},{1,0},{0,1}}]
M = S^1/(ideal {x_0^2, x_1^2, x_2^2, x_3^2})
RM = toricRR M
assert(RM.dd^2 == 0)
assert(isHomogeneous RM)
///


toricLL = method();
--Input: N a (multi)-graded E-module.
--Caveat: Assumes N is finitely generated.
--Caveat 2:  arrows of toricLL(N) correspond to exterior multiplication (not contraction)
toricLL Module := N -> (
    E := ring N;
    if not isSkewCommutative E then error "ring N is not skew commutative";
    if not E.?symmetric then E.symmetric = dualRingToric E;
    S := E.symmetric;
    N = coker presentation N;
    bb := basis N;
    b := degrees source bb;
    homDegs := sort unique apply(b, i-> last i);
    inds := new HashTable from apply(homDegs, i-> i=> select(#b, j-> last(b#j) == i));
    sBasis := new HashTable from apply(homDegs, i-> i => bb_(inds#i));
    FF := new HashTable from apply(homDegs, i ->(
	    i => S^(apply((degrees sBasis#i)_1, j-> -drop(j,-1)))
	    )
	);
    EtoS := map(S,E,toList (numgens S:0));
    differential := sum for i to numgens E-1 list (
        S_i* EtoS matrix basis(-infinity, infinity, map(N,N,E_i))
	);
    if #homDegs == 1 then (chainComplex map(S^0,FF#0,0))[1] else (
	chainComplex apply(drop(homDegs,-1), i-> map(FF#i,FF#(i+1), (-1)^(homDegs#0)*differential_(inds#(i+1))^(inds#(i))))[homDegs#0]
	)
    )

stronglyLinearStrand = method();
stronglyLinearStrand Module := M -> (
    S := ring M;
    h := heft S;
    if h === null then error("--ring M does not have heft vector");
    if not same degrees M then error("--M needs to be generated in a single degree");
    degM := first degrees M;
    canonical := sum flatten degrees vars S;
    degrange := unique prepend(degM, apply(degrees S, d -> d + degM));
    RM := toricRR(M,degrange);
    mat := RM.dd_0;
    cols := positions(degrees source mat, x -> drop(x,-1) == degM + canonical);
    N := ker mat_cols;
    toricLL N
)


TEST///
restart
load "MultigradedBGG.m2"
needsPackage "NormalToricVarieties"
S = ring hirzebruchSurface 3;
M = coker vars S;
stronglyLinearStrand(M)--should give Koszul complex, and it does.
M = coker matrix{{x_1, x_2^2}}
stronglyLinearStrand(M)--gives correct answer.
S = ring weightedProjectiveSpace {1,1,1,2,2}
I = minors(2, matrix{{x_0, x_1, x_2^2, x_3}, {x_1, x_2, x_3, x_4}})
M = Ext^3(module S/I, S^{{-7}})
N = stronglyLinearStrand M
--Wrong answer. The terms are correct, and C.dd_2 is correct.
--But C.dd_1 is, bizarrely, 0. I believe the correct matrix is:
--matrix{{-x_0, 0, -x_3, 0, x_1, 0}, {-x_1, -x_0, -x_4, -x_3, x_2, x_1}, {0, -x_1, 0, -x_4, 0, x_2}}
--I think stronglyLinearStrand is correct. The problem is with toricLL. Not sure yet what
--the problem is. In particular, I'm not sure if this is the same problem as before or a new problem.
--Aside: M is the canonical module of the coordinate ring of a copy of P^1 embedded in
--the weighted projective space P(1,1,1,2,2). 
///

TEST///
restart
load "MultigradedBGG.m2"
needsPackage "NormalToricVarieties"
S = ring hirzebruchSurface 3;
E = dualRingToric S;
C = toricLL(E^1)
C_(-1)
C.dd
C.dd
N = module ideal {e_0, e_1*e_3}

C3 = toricLL(module ideal(e_0, e_1*e_3))
C3.dd
C4 = toricLL(module ideal(e_2, e_1*e_3))
C4.dd

--silly rank 1 example
toricLL(coker vars E)
--applying LL to a rank 1 free module should give a Koszul complex (up to a degree twist)
C1 = toricLL(E^1)
isHomogeneous C1
(C1.dd)^2 == 0
C1.dd
(dual C1).dd
--how about a cyclic but non-free module:
C2 = toricLL(coker matrix{{e_0, e_1}})
isHomogeneous C2
(C2.dd)^2 == 0
C2.dd
isHomogeneous oo
--let's try it with a non-cyclic module
S = ring weightedProjectiveSpace {1,1,1,1}
E = dualRingToric S
N = module ideal(e_2, e_1*e_3)
presentation N
module ideal(e_0, e_1 * e_3)
C3 = toricLL(module ideal(e_0, e_1*e_3))
C3.dd
C4 = toricLL(module ideal(e_2, e_1*e_3))
C4.dd
isHomogeneous C3
(C3.dd)^2 == 0
C3.dd
C3_(-1)

N = coker matrix{{e_0, e_1}}
C = toricLL(N)
isHomogeneous oo
(toricLL(N)).dd
degrees C_0
degrees C_1
degrees C_2



--rank 1 example
toricLL(coker vars E)
isHomogeneous oo
--we should get the (twisted) Koszul complex when we input E
N = E^1
toricLL(E^1)
isHomogeneous oo




N' = ker matrix{{e_1}, {e_2}, {e_3}}
toricLL(N')
N'' = coker presentation N'
toricLL(N'')
///




sheafCohomologyBGG = method();
--Input: M a graded S-module, S the Cox ring of a weighted projective space.
--output will be HH^i \tilde( M(j) ). Here we're computing cohomology over the associated weighted projective stack.
sheafCohomologyBGG(Module, ZZ, ZZ) := (M, i, j) -> (
    if depth M <= 1 then error("--M is not saturated");
    if i == 0 then return M_j;
    r := regularity M;
    minDeg := min flatten degrees gens M;
    maxDeg := max flatten degrees gens M;
    --LL := toList (minDeg..maxDeg);
    LL := {{0}, {1}, {2}, {3}};
    RM := toricRR(M,LL);
    F := resDM(RM, LengthLimit => 4);
    T := F_0;
    --numRows Hom(coker vars ring T, T)_{j, -i}
    basis({-6,-3}, Hom(coker vars ring T, T));
    )



TEST///
restart
needsPackage "NormalToricVarieties"
load "MultigradedBGG.m2"
X = weightedProjectiveSpace {1,1,2}
S = ring X
M = S^1
regularity M
M = S^{{6}}
regularity M
M = coker vars S
regularity M
needsPackage "Depth"
sheafCohomologyBGG(S^1,1,-6)
///

end;

--TESTS


--DEMO. We should delete this eventually, but leaving it here for now in case it's useful.
restart
load "MultigradedBGG.m2"
--We introduce a new type: a differential module.
--There is a method that turns any map A such that A^2 = 0 into a differential module.
--Here is an example from the talk over the weekend:
R = ZZ/101[x, y]
M = map(R^2, R^2, matrix{{x*y, -x^2}, {y^2, -x*y}})
D = differentialModule(M)
oo.dd
--We can compute its homology, as in the exercises:
mingens HH(D)
--We can take its free resolution in two ways. First, we execute
--the algorithm from the problem set over the weekend ("KC" stands for "killing cycles")
resKC(D)
oo.dd
--Here is another method for resolving. The underlying module of the output will always be
--the free resolution of H(D):
resDM(D)
oo.dd
--(it gives the same answer as resKC in this example, but this need not always be the case)
--This resolution isn't minimal. Let's minimize it:
minimizeDM resDM(D)
oo.dd
--We get back what we started, as expected.

--Application to BGG:
--The output of the multigraded BGG functor is a differential module.
--Here is an example from the problem set over the weekend:
restart
load "MultigradedBGG.m2"
X = weightedProjectiveSpace {1,1,2}
S = ring X
--There is a BGG functor R: mod(S) --> DM(E),
--where E is the "Koszul dual" of S. It induces an equivalence on derived categories.
--The output of R is typically an infinitely generated free E-module, so we must
--specify a degree range for the output:
toricRR(S^1, {0, 1, 2, 3})
oo.dd
--We also implement the left adjoint, L, of R. Applying it to a rank 1 free module gives
--the Koszul complex:
E = dualRingToric S
toricLL(E^1)
oo.dd
--Given an S-module M, we can compute H^i(X,\widetilde{M}(j)) using BGG. Here is an example.
M = S^1
--Compute R(M) (in a certain degree range):
 = toricRR(M, {0,1,2,3})
--Resolve it:
tail = (resDM(D, LengthLimit => 3))
--each H^2(X, O(j)) can be computed by picking off the socle generators of the
--summands of the tail
bggCohomologyCheck j -> (
    kk = coker vars ring tail_0;
    sum flatten entries basis({j,-3}, Hom(kk, tail_0)) == rank HH^2(X, sheaf(S^{{j}}))
    )    
--bggCohomology prints "true" if the rank of the degree (j, -3) part of 
--the socle of the free resolution of R(M) is equal to dim_k H^2(X, O(j))
for k from -7 to 7 do (
    print bggCohomologyCheck(k);
    );
--we get the right answer in a certain "window".
bggCohomologyCheck(-8)
--we get the wrong answer eventually, because our window for R(M) and its resolution 
--is only so big. To get more cohomology, one must widen the window.

