-*
newPackage("MultigradedBGG",
    Version => "1.1",
    Date => "5 June 2023",
    Headline => "Computing with the multigraded BGG correspondence",
    Authors => {
        {Name => "Michael K. Brown",         Email => "mkb0096@auburn.edu",    HomePage => "http://webhome.auburn.edu/~mkb0096/" }
	{Name => "Daniel Erman",    	     Email => "erman@wisc.edu",        HomePage => "https://people.math.wisc.edu/~erman/" }
	{Name => "Tara Gomes",	    	     Email => "gomes072@umn.edu",      HomePage => "Fill in" }
	{Name => "Pouya	Layeghi",    	     Email => "layeg001@umn.edu",      HomePage => "Fill in" }
	{Name => "Prashanth Sridhar",	     Email => "pzs0094@auburn.ed",     HomePage => "https://sites.google.com/view/prashanthsridhar/home" }
	{Name => "Andrew Tawfeek",   	     Email => "atawfeek@uw.edu",       HomePage => "https://www.atawfeek.com/" }
	{Name => "Eduardo Torres Davila",    Email => "torre680@umn.edu",      HomePage => "https://etdavila10.github.io/" }
	{Name => "Jay Yang",         	     Email => "jayy@wustl.edu",        HomePage => "https://www.math.wustl.edu/~jayy/" }
	{Name => "Sasha	Zotine",    	     Email => "18az45@qiueensu.ca",    HomePage => "Fill in" }
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


needsPackage "Polyhedra";
load "DifferentialModules.m2";
---
---
---
---

--Input:    A pair of matrices (M,N)
--Output:   The effect of contracting M by N. 
matrixContract=method()
matrixContract(Matrix,Matrix) := (M,N) -> (
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

dualRingToric(PolynomialRing) := opts -> (S) ->(
--  Input: a polynomial ring OR an exterior algebra
--  Output:  the Koszul dual exterior algebra, but with an additional 
--           ZZ-degree, the ``standard grading'' where elements of \bigwedge^k
--           have degree k.
    kk := coefficientRing S;
    degs := null;
    if isSkewCommutative S == false then(
    	degs = apply(degrees S,d-> (-d)|{1});
	e := opts.SkewVariable;
    	ee := apply(#gens S, i-> e_i);
    	return kk[ee,Degrees=>degs,SkewCommutative=>true]
	); 
    if isSkewCommutative S == true then(
    	degs = apply(degrees S,d-> remove(-d,-1));
	y := opts.Variable;
    	yy := apply(#gens S, i-> y_i);
    	return kk[yy,Degrees=>degs,SkewCommutative=>false]
	);       
    );

TEST ///
restart
needsPackage "NormalToricVarieties"
S = ring(hirzebruchSurface(2, Variable => y));
E = dualRingToric(S,SkewVariable => f);
SY = dualRingToric(E);
assert(degrees SY == degrees S)
///

--I think toricRR works well
toricRR = method();
--Input: (M,LL) M a (multi)-graded S-module.
--             LL a list of degrees.
--Output:  The differenial module RR(M) in degrees from LL,
--         presented as a complex in homological degrees -1, 0 ,1
--         and with the same differential in both spots.
toricRR(Module,List) := (M,LL) ->(
    S := ring(M);
    if not isCommutative S then error "ring M is not commutative";
    if not S.?exterior then S.exterior = dualRingToric(S);
    E := S.exterior;
    relationsM := presentation M;
    -- this used to say "gens image presentation M"... just in case a bug arises
    f0 := gens image basis(LL_0,M);
    scan(#LL-1,i-> f0 = f0 | gens image basis(LL_(i+1),M));
    df0 := apply(degrees source f0,i-> (-1)*i|{0});
    df1 := apply(degrees source f0,i-> (-1)*i|{-1});
    dfneg1 := apply(degrees source f0,i-> (-1)*i|{1});
    SE := S**E;
    --the line below is better for degrees,it overwrites S somehow...
    --SE := coefficientRing(S)[gens S|gens E, Degrees => apply(degrees S,d->d|{0}) | degrees E, SkewCommutative => gens E];
    tr := sum(dim S, i-> SE_i*SE_(dim S+i));
    newf0 := sub(f0,SE)*tr;
    relationsMinSE := sub(relationsM,SE);
    newf0 = newf0 % relationsMinSE;
    newg := matrixContract(transpose sub(f0,SE),newf0);
    g' := sub(newg,E);
    differentialModule(chainComplex{map(E^dfneg1,E^df0, g'),map(E^df0,E^df1, g')}[1])
    )



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

--I think toricLL doesn't work and needs to be debugged.
toricLL = method();
--Input: N a (multi)-graded E-module.
--Caveat: Assumes N is finitely generated.
--Caveat 2:  arrows of toricLL(N) correspond to exterior multiplication (not contraction)
toricLL(Module) := (N) ->(
    E := ring(N);
    if not isSkewCommutative E then error "ring N is not skew commutative";
    if not E.?symmetric then E.symmetric = dualRingToric(E);
    S := E.symmetric;
    bb := basis(N);
    b := (degrees source bb);
    homDegs := sort unique apply(b, i-> last i);
    inds := new HashTable from apply(homDegs, i-> i=> select(#b, j-> last(b#j) == i));
    sBasis := new HashTable from apply(homDegs, i-> i => (bb)_(inds#i));
    FF := new HashTable from apply(homDegs, i ->(
	    i => S^(apply((degrees sBasis#i)_1, j-> remove(j,-1)))
	    )
	);
    relationsN := presentation N;
    SE := S**E;
    --SE := coefficientRing(S)[gens S|gens E, Degrees => apply(degrees S,d->d|{0}) | degrees E, SkewCommutative => gens E];
    --why is this overwriting the definition of e_i?
    tr := sum(dim S, i-> SE_i*SE_(dim S+i));
    f0 := gens image basis(N);
    newf0 := sub(f0,SE)*tr;
    relationsMinSE := sub(relationsN,SE);
    newf0 = newf0 % relationsMinSE;
    newg := matrixContract(transpose sub(f0,SE),newf0);
    g' := sub(newg,S);
    --Now we have to pick up pieces of g' and put them in the right homological degree.
    --Note:  perhaps we want everything transposed??
    dual(chainComplex apply(remove(homDegs,-1), i-> map(FF#i,FF#(i+1),transpose g'_(inds#i)^(inds#(i+1))))[-homDegs#0])
    )


-* The following commented-out code may be useful, but I can't tell if we'll need it yet.
--I think it was originally intended for applications that we're not going to pursue.


--  Input:  a polynomial ring S
--  Output:  a list of all degrees between 0 and the anti-canonical degree
--  CAVEAT:  needs "Polyhedra" package.
koszulDegrees = method();
koszulDegrees(PolynomialRing) := S ->(
    K := koszul vars S;
    L := unique flatten apply(dim S+1,i-> degrees K_i);
    P := convexHull transpose matrix L;
    LP := latticePoints P;
    sort flatten apply(#LP,i-> entries transpose (latticePoints P)_(i))
    )

--  Input:  a polynomial ring S
--  Output:  a list of all degrees between 0 and the anti-canonical degree
--  CAVEAT:  needs "Polyhedra" package.
effectivePart = method();
effectivePart(PolynomialRing,ZZ) := (S,n) ->(
    zd := {apply(#(degree S_0),i-> 0)};
    ed := degrees S|zd;
    P := convexHull transpose (n*(matrix ed));
    LP := latticePoints P;
    sort flatten apply(#LP,i-> entries transpose (latticePoints P)_(i))
    )

--  Input:  degrees d1, d2, and a ring S
--          Assumes d1,d2 are lists of integers of length n and that the variables
--          of S have degrees which are lists of integers of length n.
--  Output:  Checks whether d1 - d2 is an effective degree on S.
greaterEqual=method()
greaterEqual(List,List,Ring) := (d1,d2,S) -> hilbertFunction(d1-d2,S) != 0


--  Input:  A module/matrix/chain complex and a multidegree d.
--  Output:  The truncated module/matrix/chain complex in degrees <= d.
degreeTruncation=method()
degreeTruncation(Matrix,List) := (M,d) -> (
    --rows and cols with degrees <= d.
    S := ring(M);
    rows:= positions(degrees target M,d'-> greaterEqual(d,d',S));
    columns:= positions(degrees source M,d'-> greaterEqual(d,d',S));
    sourceInc := (id_(source M))_columns;
    targetInc := (id_(target M))_rows;    
    ((M^rows)_columns,targetInc,sourceInc)
    )

degreeTruncation(ChainComplex,List) := (K,d) -> (
    --subcomplex where all rows and cols of differentials have degrees <= d.
    S := ring(K);
    a := min K;
    Ka := K[a];
    L := apply(length Ka,i->degreeTruncation(Ka.dd_(i+1),d));
--    L := apply(toList(min(K[a]).. max(K[a])-1),i->degreeTruncation(Ka.dd_(i+1),d));    
    tKa:=chainComplex(apply(length Ka,i->L_i_0));
    phi := map(K,tKa[-a], i-> if i == a then L_(i-a)_1 else L_(i-a-1)_2);
--    assert(isChainComplexMap phi);
    phi
    )


degreeTruncation(ChainComplexMap,List) := (phi,d) -> (
    S := ring(phi);
    G := target phi;
    F := source phi;
    phiF := degreeTruncation(F,d);
    F' := source (phiF = degreeTruncation(F,d));
    G' := source (phiG = degreeTruncation(G,d));
    map(G', F', i->(phi_i * phiF_i)// phiG_i)
    )


--   Input:  a ring S and a multidegree d.
--   Output:  the truncated Koszul complex of summands of degrees <= d.
truncatedKoszul = method();
truncatedKoszul(Ring,List) := (S,d)->(
    E := dualRingToric S;
    F := (toricLL(E^1))[-dim S];
    canon := sum degrees S;
    G := F**(ring F)^{-canon+d};
    source degreeTruncation(G,{0})
    )


--   Input:  a ring S, degrees d1 and d2, and a map f between koszul complexes of degree d2-d1.
--   Output:  the inducted map on truncated Koszul complexes.
truncatedKoszulMap = method();
truncatedKoszulMap(Ring,List,List,RingElement) := (S,d1,d2,f)->(
    if -remove(degree f,-1) != d2-d1 then error "Wrong degree for such a map";
    N1 := ptrunc(E^1,d1);
    N2 := ptrunc(E^1,d2);
    phi := map(N1**E^{d1|{0}},N2**E^{d2|{-last degree f}},matrix{{f}});
    toricLL phi
    )

*-


--This next bit is for the sheaf cohomology algorithm. It gives code for resolving RR(M).
--We will only use it in the case of a Z-grading, which is far simpler than what is being
--attempted here. Most of this can probably be removed, but I left it here in case it's
--helpful.

--What follows are old comments:
--The following develops a "corner complex" code.
--For weighted projective spaces, this can be used to compute parts of the Tate resolution,
--and therefore compute sheaf cohomology.
--But for other toric varieties, these corner complexes seem to compute something else.
--This actually raises a significant theoretical question:  is there any analogue of the corner complex
--for arbitrary toric varieties??

--Input:  a differential module F, and a corner degree cDeg.
--Output:  a new differential module F' where we have adjoined "killed" cycles of F of degree <= cDeg
oneStepCorner = method();
oneStepCorner(List,DifferentialModule) := (cDeg,F) -> (
    E := ring F;
    dF := degree F;
    H := res(coker F.dd_1,LengthLimit => 2);
    newSyz := (H.dd_2**(ring H)^{dF}) % H.dd_1;
    G = chainComplex newSyz;
    newSpots = {};
    D := degrees G_1;
    D1 := apply(D,i-> remove(i,-1));
--  keep nontrivial new syzygies which are "below" the corner
    scan(#D1, i->(if greaterEqual(cDeg,D1#i,S) and D1#i != cDeg and submatrix(newSyz,{i}) != 0 
	    then(newSpots = newSpots|{i};)));
    Gpre := E^(apply(newSpots,i-> -D#i + dF));
    phi := map(G_0,Gpre,(newSyz)_(newSpots));
    psi := gens trim image(phi % (matrix H.dd_1));
    Gnew := source psi;
    d1 := psi|H.dd_1;    
    d2 := map(Gnew**E^{dF},Gnew++H_1, (i,j)->0);
    d := (d2||d1);
    differentialModule(chainComplex(d**E^{dF},d)[1])
)




cornerDM = method(TypicalValue => DifferentialModule,
            Options => {
    	    LengthLimit => 3
    	    });


--Input:  a differential module F, and a corner degree cDeg, and an (optional) LengthLimit.
--Output:  a new differential module F' where we have adjoined "killed" cycles of F of degree <= cDeg,
--         and we have itereated that procedure for the specified number of steps.
cornerDM(List,DifferentialModule) := opts -> (cDeg,F)->(
    n := opts.LengthLimit;
    scan(n, i-> F = oneStepCorner(cDeg,F));
    F
    )   



--
--Specialized functions for a Hirzebruch surface.
--


---
--  Does oneStepHirz but also lifts the map phi to the new module.
oneStepRows = (topDeg,F,rows,aut) -> (
    E := ring F;
    dF := degree F;
    H := res(coker F.dd_1,LengthLimit => 2);
    newSyz := (H.dd_2**(ring H)^{dF}) % H.dd_1;
    G = chainComplex newSyz;
    newSpots = {};
    D := degrees G_1;
    D1 := apply(D,i-> remove(i,-1));
--  We keep the syzygies which are:
--    "below" a specified degree topDeg (to avoid boundary effects)
--    in the specified rows, and
--    are nontrivial.
    scan(#D1, i->(if greaterEqual(topDeg,D1#i,S) and isSubset(set{D1#i#1},set(rows)) and submatrix(newSyz,{i}) != 0 
	    then(newSpots = newSpots|{i};)));
    phi := map(G_0,E^(apply(newSpots,i-> -D#i + dF)),(newSyz)_(newSpots));
--    assert(H.dd_1*phi == 0);
    psi := gens trim image(phi % (matrix H.dd_1));
--    assert(H.dd_1*psi == 0);
    Gnew := source psi;
    psiGF := map(Gnew**E^{dF},Gnew,0) || psi;
--This maps the new syzygies into G++F where G=Gnew and F=F_0.
    d1 := psi|(-1)*H.dd_1; 
    d2 := map(Gnew**E^{dF},Gnew++H_1, (i,j)->0);
    d := (d2||d1);
--This is the new differential on G++F.
--    assert(d^2 == 0)
--Now we extend the automorphism.
    oldAutGF := map(Gnew**E^{dF}, Gnew, (i,j) ->0)++aut;
    newAutGF := map(source d**E^{dF},Gnew,aut*psi // d1);
    --newAutGF := map(source d**E^{dF},Gnew,oldAutGF*(psiGF**E^{-dF}) // d);
    newAut := newAutGF | (map(Gnew,source aut, (i,j)->0)||-aut);
--    assert(-d*newAut - newAut*d == 0)
--assert(newAut^2 == 0).. 
    newF := differentialModule(chainComplex(d**E^{dF},d)[1]);
    (newF, newAut)  
)

multiStepRows = (topDeg,F,rows,aut,n) -> (
    for i from 1 to n do(
	LL := oneStepRows(topDeg,F,rows,aut);
	F = LL_0;
	aut = LL_1;
	);
    (F,aut)
)

--assumes "rows" is a singleton set.
oneStepVert = (botLeft,topRight,F,rows) -> (
    E := ring F;
--This code uses the parameter of the hirzebruch.
    a := (degree E_1)_0;
    dF := degree F;
    H := res(coker F.dd_1,LengthLimit => 2);
    newSyz := (H.dd_2**(ring H)^{dF}) % H.dd_1;
    G = chainComplex newSyz;
    newSpots = {};
    D := degrees G_1;
    D1 := apply(D,i-> remove(i,-1));
--    We keep the syzygies which are:
--    in the triangle (directions spanned by degrees of e_1,e_3)
--    in the specified rows, and are nontrivial.
--    scan(#D1, i->(if D1#i#0+D1#i#1 >= sum botLeft  and isSubset(set{D1#i#1},set(rows)) and(
    scan(#D1, i->(if D1#i#0+a*(D1#i#1) >= sum botLeft  and D1#i#1 <= min rows+1 and(
		 submatrix(newSyz,{i}) != 0 and D1#i#0 <= topRight#0 )
	    then(newSpots = newSpots|{i};)));
    phi := map(G_0,E^(apply(newSpots,i-> -D#i + dF)),(newSyz)_(newSpots));
--    assert(H.dd_1*phi == 0);
    psi := gens trim image(phi % (matrix H.dd_1));
--    assert(H.dd_1*psi == 0);
    Gnew := source psi;
    psiGF := map(Gnew**E^{dF},Gnew,0) || psi;
--This maps the new syzygies into G++F where G=Gnew and F=F_0.
    d1 := psi|(-1)*H.dd_1; 
    d2 := map(Gnew**E^{dF},Gnew++H_1, (i,j)->0);
    d := (d2||d1);
    differentialModule(chainComplex(d**E^{dF},d)[1])
    )


multiStepVert = (botLeft,topRight,F,rows,n)->(
     for i from 1 to n do(
	 F = oneStepVert(botLeft,topRight,F,rows);
	 --rows = rows + {-1};
	 );
    F
    )

botLeftAndTopRight = L->(
    min1 = min apply(L,l-> l_1);
    bL = min select(L,l-> l_1 == min1);
    max1 = max apply(L,l-> l_1);
    tR = max select(L,l-> l_1 == max1);
    (bL,tR)
    )

topLeft = L->(
    max1 = max apply(L,l-> l_1);
    tL = min select(L,l-> l_1 == max1);
    tL
    )

hirzResolve = (M,L,n1,n2)->(
    rows := unique apply(L,i-> i_1);
    F := toricRR(M,L,{0,2});
    aut = differential toricRR(M,L,{1,3});
    topDeg := topLeft(L);
    GA := multiStepRows(topDeg,F,rows,aut,n1);
    vertF := differentialModule( chainComplex(GA_1**(ring (GA_0))^{degree (GA_0)},GA_1)[1]);
    (botLeft,topRight) = botLeftAndTopRight(unique degrees vertF_0);
    multiStepVert(botLeft,topRight,vertF,{min rows - 1},n2)
    )

--Input:  the list of degrees computed in a Tate resolution.
--Output:  the vertices of the corresponding polytope.
corners = L->(
    LL := apply(L,l-> remove(l,-1));
    P := convexHull transpose matrix LL;
    vertices P
    )

matrixCorners = L->(
    LL := apply(L,l-> remove(l,-1));
    min1 := min apply(L,l-> l_1);
    max1 := max apply(L,l-> l_1);
    min0 := min apply(L,l-> l_0);
    max0 := max apply(L,l-> l_0);
    a := max0 - min0 + 1;
    b := max1 - min1 + 1;
    transpose map(ZZ^a,ZZ^b, (i,j) -> if (i+min0,max1 - j) == (0,0) then 5 else(
	    if isSubset(set({{i+min0,max1 - j}}),set(LL)) then 1 else 0)
	)
    )


makeConvex = L->(
    needsPackage "Polyhedra";
    P = convexHull transpose matrix L;
    LP = latticePoints P;
    flatten apply(#LP,i-> entries transpose (latticePoints P)_(i))
    )

end;


--TESTS
restart
load "MultigradedBGG.m2"
loadPackage "NormalToricVarieties"
X = weightedProjectiveSpace {1,1,1}
S = ring X
E = dualRingToric S
toricLL(E^1)
N = coker vars E
toricLL(N)









