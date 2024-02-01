-*
newPackage("MultigradedBGG",
    Version => "1.1",
    Date => "5 June 2023",
    Headline => "The multigraded BGG correspondence in Macaulay2",
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


exports {
    "dualRingToric",
    "toricRR"
    "toricLL"
    "stronglyLinearStrand"
    "unfold"
    "DifferentialModule",
    "differentialModule",
    "unfold",
    "foldComplex",
    "resDM",  
    "resKC",    
    "minimizeDM",
    "differential"
    }
*-


--path=prepend("../",path)
--Input:    A list of matrices with the same number of rows.
--Output:   The concatenation of those matrices.
concatMatrices=method()
concatMatrices(List) := L -> (
    m:= first L;
    scan(#L-1,i->m=m|L_(i+1));
    m)


--Example
--R = ZZ/101[x,y]
--A = matrix{{x*y, -x^2}, {y^2, -x*y}}
--d = map(R^2, R^2, A)
--C = chainComplex(d, d)


---
---
---
---GENERAL STUFF ON DIFFERENTIAL MODULES
---
---
---
DifferentialModule = new Type of ChainComplex
DifferentialModule.synonym = "differential module"


differentialModule = method(TypicalValue => DifferentialModule)
differentialModule ChainComplex := C -> (
    --add error if C is not of the following form:
    --C_1 --> C_0 --> C_{-1}, where the two differentials
    --are the same matrix (say of degree d, if there 
    --is a grading), this matrix squares to 0, C_1 = C_0(-d), and C_{-1} = C_0(d).
  new DifferentialModule from C);


---MAYA: changed this so that source and target are the same, map may be nonzero degree
differentialModule Matrix := phi -> (
    --check if the source and target are the same up to a twist
    if phi^2 != 0 then error "The differential does not square to zero.";
    R := ring phi;
    -- MAYA d := (degrees source phi)_0 - (degrees target phi)_0;
    -- MAYA if target phi != source phi**R^{d} then error "source and target of map are not the same, up to a twist"; 
    if target phi != source phi then error "source and target of map are not the same"; 
    -- MAYA new DifferentialModule from (chainComplex(phi**R^{d},phi)[1]));
    new DifferentialModule from (chainComplex(phi,phi)[1]));

-*
differentialModule Matrix := phi -> (
    --check if the source and target are the same up to a twist
    R := ring phi;
    d := (degrees source phi)_0 - (degrees target phi)_0;
    if target phi != source phi**R^{d} then error "source and target of map are not the same, up to a twist"; 
    --if target phi != source phi then error "source and target of map are not the same"; 
    new DifferentialModule from (chainComplex(phi**R^{d},phi)[1]));
    --new DifferentialModule from (chainComplex(phi,phi)[1]));
 *-


ring(DifferentialModule) := Ring => D -> D.ring;
module DifferentialModule :=  (cacheValue symbol module)(D -> D_0);
--MAYA degree DifferentialModule := ZZ => (D -> (degrees D_1)_0 - (degrees D_0)_0);
degree DifferentialModule := ZZ => (D -> degree D.dd_1); --maya changed to make degree of dm the degree of the map
differential = method();
differential DifferentialModule := Matrix=> (D->D.dd_1);
kernel DifferentialModule := Module => opts -> (D -> kernel D.dd_0); 
image DifferentialModule := Module => (D -> image D.dd_1); 
homology DifferentialModule := Module => opts -> (D -> HH_0 D);

isFreeModule(DifferentialModule) := D ->(
    isFreeModule module D
    )

--MAYA: changed to be compatible with new degree convention
unfold = method();
--Input:  a differential module and a pair of integers low and high
--Output:  the unfolded chain complex of the differential module, in homological degrees
--         low through high.
unfold(DifferentialModule,ZZ,ZZ) := ChainComplex => (D,low,high)->(
    L := toList(low..high);
    d := degree D;
    R := ring D;
    phi := differential D;
    --MAYA chainComplex apply(L,l-> phi**R^{l*d})[-low]
    chainComplex apply(L,l-> phi)[-low]
    )


minFlagOneStep = method()
minFlagOneStep(DifferentialModule) := (D) -> (
    d := degree D;
    R := ring D;
    minDegree := min degrees trim HH_0 D;
    colList := select(rank source (mingens HH_0(D)), i -> degree (mingens HH_0(D))_i == minDegree);
    minDegHom := (mingens HH_0(D))_colList;
    homMat :=  mingens image((minDegHom) % (image D.dd_1));
    G := res image homMat;
    psi := map(D_0,G_0**R^{d},homMat, Degree=>d);
    newDiff := matrix{{D.dd_1,psi},{map(G_0**R^{d},D_1,0, Degree=>d), map(G_0**R^{d},G_0**R^{d},0, Degree=>d)}};
    assert (newDiff*(newDiff) == 0);
    differentialModule newDiff
)

resMinFlag = method();
resMinFlag(DifferentialModule, ZZ) := (D, k) -> (
    d := degree D;
    R := ring D;
    s := numgens D_0;
    scan(k,i-> D = minFlagOneStep(D));
    t := numgens D_0;
    newDiff := submatrix(D.dd_1, toList(s..t-1),toList(s..t-1));
    differentialModule (chainComplex(newDiff**R^{-d},newDiff**R^{-d})[1])
)

-- MAYA: changed to be compatible with new degree convention         
killingCyclesOneStep = method();
killingCyclesOneStep(DifferentialModule) := (D)->(
--  commented out code was working but I'm trying to improve
--    G := res HH_0 D;
--    psi := map(D_0,G_0,(gens (HH_0 D)) %  (image D.dd_1));
--    d := degree D;
--    R := ring D;
--    newDiff := matrix{{D.dd_1,psi},{map(G_0**R^{d},D_1,0),map(G_0**R^{d},G_0,0)}};
--    assert (newDiff*(newDiff**R^{-d}) == 0);
--    differentialModule newDiff
--    )
    d := degree D;
    R := ring D;
    homMat := mingens image((gens HH_0 D) %  (image D.dd_1));
    G := res image homMat;
    psi := map(D_0,G_0**R^{d},homMat, Degree=>d);
    --MAYA newDiff := matrix{{D.dd_1,psi},{map(G_0**R^{d},D_1,0),map(G_0**R^{d},G_0,0)}};
    newDiff := matrix{{D.dd_1,psi},{map(G_0**R^{d},D_1,0, Degree=>d),map(G_0**R^{d},G_0**R^{d},0, Degree=>d)}}; --maya added
    assert (newDiff*(newDiff) == 0);
    differentialModule newDiff
    )
--killing cycles resolution
--Input:  a free(?) differential module D and an integer k
--Output:  the killing cycles resolution of D after k steps.
-- NOTE:  until recently this produced the CONE over the resolution.
--        but I'm interested in the resolution itself so the last few lines change that
--MAYA: changed to be compatible with new degree convention
--Michael: changed the name of this method (and the next) to resDM, rather than resKC. The old
--resDM has been removed.

resDM = method();
resDM(DifferentialModule,ZZ) := (D,k)->(
    d := degree D;
    R := ring D;
    s := numgens D_0;
    scan(k,i-> D = killingCyclesOneStep(D));
    t := numgens D_0;
    newDiff := submatrix(D.dd_1, toList(s..t-1),toList(s..t-1));
    differentialModule (chainComplex(newDiff**R^{-d},newDiff**R^{-d})[1])
    )

--same as above, but default value of k 
--MAYA: changed to be compatible with new degree convention
resDM(DifferentialModule) := (D)->(
    k := dim ring D + 1;
    d := degree D;
    R := ring D;
    s := numgens D_0;
    scan(k,i-> D = killingCyclesOneStep(D));
    t := numgens D_0;
    newDiff := submatrix(D.dd_1, toList(s..t-1),toList(s..t-1));
    differentialModule (chainComplex(newDiff**R^{-d},newDiff**R^{-d})[1])
    )



--  Subroutines and routines to produce the minimal part of a matrix.
--  NEEDS TO BE A SQUARE MATRIX
--  MAYA: matrix should be a map of free modules
-- MAYA: changed to preserve degree
minimizeDiffOnce = method();
minimizeDiffOnce(Matrix,ZZ,ZZ) := (A,u,v) -> (
    a := rank target A;
    R := ring A;
    inv := (A_(u,v))^(-1);
    N := map(source A,source A, (i,j) -> if i == v and j != v then -inv*A_(u,j) else 0) + id_(R^a);
    Q := map(target A,target A, (i,j) -> if j == u and i != u then -inv*A_(i,v) else 0) + id_(R^a);
    A' := Q*N^(-1)*A*N*Q^(-1);
    newRows := select(a,i-> i != u and i != v);
    newCols := select(a,i-> i != u and i != v);
    A'_(newCols)^(newRows)
    )

units = method();
units(Matrix) := A->(
    a := rank source A;
    L := select((0,0)..(a-1,a-1), (i,j)->  isUnit A_(i,j));
    L
    )

--  Input: A SQUARE matrix
--  Output: a minimization of it.
minimizeDiff = method();
minimizeDiff(Matrix) := A ->(
    NU := #(units A);
    while NU > 0 do(
 	L := units A;
	(u,v) := L_0;
	A = minimizeDiffOnce(A,u,v);
	NU = #(units A);
	);
    A
    )


--  Input:  A (finite, free) differential module
--  Output: A minimization of that DM.
-- MAYA: changed to be compatible with new degree convention
minimizeDM = method();
minimizeDM(DifferentialModule) := r ->(
    R := ring r;
    d := degree r;
    A := minimizeDiff(r.dd_1);
    degA := map(target A, source A, A, Degree=>d);
    differentialModule (chainComplex(degA,degA)[1])
    )

---
---

--  Input:  a free complex F and a degree d
--  Output: the corresponding free differential module of degree da
--MAYA: changed to be compatible with degree convention
foldComplex = method();
foldComplex(ChainComplex,ZZ) := DifferentialModule => (F,d)->(
    R := ring F;
    L := apply(length F+1,j->(
	    transpose concatMatrices apply(length F+1,i->map(F_j,F_i, if i == j+1 then F.dd_i else 0))
	));
    FDiff := transpose(concatMatrices L);
    FMod := F_0;
    scan(length F+1, i-> FMod = FMod ++ ((F_(i+1))**(R)^{(i+1)*d}));
    degFDiff := map(FMod,FMod,FDiff, Degree=>d); --maya added
    -- MAYA differentialModule(chainComplex(FDiff**(ring F)^{d},FDiff)[1])
    differentialModule(chainComplex(degFDiff,degFDiff)[1]) --maya added
    )

-*
beginDocumentation()

doc ///
   Key 
      DifferentialModules
   Headline 
      Package for Computing Free Resolutions of Differential Modules
   Description
    Text
      Blah blah blah      
///


doc ///
   Key 
    differentialModule
    (differentialModule, map)
   Headline
    converts a square zero matrix into a differential module
   Usage
    differentialModule(f)
   Inputs
    f: module map with the same source and target
   Outputs
    : DifferentialModule 
   Description
    Text
      Given a module $f: M\to M$ of degree a this creates a degree a differential module from
      f represented as as 3-term chain complex in degree -1, 0, 1. If you want a nonzero, 
      you should specify the degree of the map explicitly.
      An error is returned if the source and target of f are not equal.
    Example
      R = QQ[x]
      phi = map(R^1/(x^2),R^1/(x^2),x, Degree=>1)
      differentialModule(phi)
///


doc ///
   Key 
    unfold
    (unfold,DifferentialModule,ZZ,ZZ)
   Headline
    converts a differential module into a 1-periodic complex
   Usage
    unfold(D,a,b)
   Inputs
    D: DifferentialModule
    a: ZZ
    b: ZZ
   Outputs
    : ChainComplex
   Description
    Text
      Given a differential module D and an integers a and b it produces a
      chain complex with the module D in homological a through b and where
      all maps are the differential of D.
    Example
      phi = matrix{{0,1},{0,0}};
      D = differentialModule(phi);
      unfold(D,-3,4)
///


doc ///
   Key 
    foldComplex
    (foldComplex,ChainComplex,ZZ)
   Headline
    converts a chain complex into a differential module
   Usage
    foldComplex(C)
   Inputs
    C: ChainComplex
    d: ZZ
   Outputs
    : DifferentialModule
   Description
    Text
      Given a chain complex C and integer d it creates the corresponding
      (flag) differential module of degree d.
    Example
      R = QQ[x,y];
      C = res ideal(x,y)
      D = foldComplex(C,0);
      D.dd_1
///


doc ///
   Key 
    resDM
    (resDM,DifferentialModule)
   Headline
    uses the Cartan-Eilenberg style construction to find a free resolution of a differential module
   Usage
    resDM(D)
   Inputs
    D: DifferentialModule
   Outputs
    : DifferentialModule
   Description
    Text
      Given a differential module D it creates a free flag resolution of D, using a Cartan-Eilenberg
      construction, up to the optional LengthLimit. So you should double
      check that.
    Example
      R = QQ[x,y];
      M = R^1/ideal(x^2,y^2);
      phi = map(M,M,x*y, Degree=>2);
      D = differentialModule phi;
      r = resDM(D)
      r.dd_1
    Text
      The default LengthLimit is 3 because I couldn't figure out how to change it. So if your 
      ring has dimension greater than 3, or if it's not a regular, then you can increase the
      LengthLimit to get more information.
    Example
      R = QQ[x]/(x^3);
      phi = map(R^1,R^1,x^2,Degree=>2);
      D = differentialModule phi;
      r = resDM(D)
      r.dd_1      
      r = resDM(D,LengthLimit => 6)
      r.dd_1      
///


doc ///
   Key 
    resKC
    (resKC,DifferentialModule)
    (resKC,DifferentialModule,ZZ)
   Headline
    uses a killing cycles style construction to find a free resolution of a differential module
   Usage
    resKC(D)
    resKC(D,k)
   Inputs
    D: DifferentialModule
    k: ZZ
   Outputs
    : DifferentialModule
   Description
    Text
      Given a differential module D it creates a free flag resolution of D, using a killing cycles
      construction.  Because of issues with adding options, there are two chocies.  The defualt
      resKC(D) runs the algorithm for the number of steps determined by the dimension of the ambient ring.
      resKC(D,k) for k steps.
    Example
      R = QQ[x,y];
      M = R^1/ideal(x^2,y^2);
      phi = map(M,M,x*y, Degree=>2);
      D = differentialModule phi;
      r = resKC(D)
      r.dd_1
    Text
      The default algorithm runs for dim R + 1 steps, again because I couldn't figure out the options.
      Adding the number of steps as a second argument is like adding a LengthLimit.
    Example
      R = QQ[x]/(x^3);
      phi = map(R^1,R^1,x^2, Degree=>2);
      D = differentialModule phi;
      r = resKC(D)
      r.dd_1      
      r = resKC(D,6)
      r.dd_1      
///


doc ///
   Key 
    DifferentialModule
   Headline
    The class of differential modules.
   Description
    Text
      A differential module is just a module with a square zero endomorphism.  We represent this
      via a 3-term complex with the module in positions -1, 0, and 1 and both maps being the
      endomorphism.
///


doc ///
   Key 
    minimizeDM
    (minimizeDM,DifferentialModule)
   Headline
    minimizes a sqaure matrix or a differential module
   Usage
    minimizeDM(D)
   Inputs
    D: DifferentialModule
   Outputs
    : DifferentialModule
   Description
    Text
      Given a differential module D this code breaks off trivial
      blocks, producing a quasi-isomorphic differential module D' with a minimal
      differential.
    Example
      R = QQ[x,y];
      M = R^1/ideal(x^2,y^2);
      phi = map(M,M,x*y, Degree=>2);
      D = differentialModule phi;
      r = resKC(D)
      r.dd_1
      mr = minimizeDM(r)
      mr.dd_1   
///


TEST /// --test basic diff mod stuff 
    S = QQ[x,y]
    m = matrix{{0,x,y,1},{0,0,0,-y},{0,0,0,x},{0,0,0,0}}
    phi = map(S^{0,1,1,2}, S^{0,1,1,2} ,m, Degree=>2)
    D = differentialModule phi
    assert(D.dd_0^2==0)
    assert(isHomogeneous D.dd_0)
    assert(degree D=={2})
    assert(prune homology D==cokernel matrix{{x,y}})
///

TEST /// --test basic diff mod stuff 2
    S = QQ[x,y]
    m = matrix{{0,x^2,x*y,1},{0,0,0,-y},{0,0,0,x},{0,0,0,0}}
    phi = map(S^{0,1,1,3}, S^{0,1,1,3} ,m, Degree=>3)
    D = differentialModule phi
    assert(D.dd_0^2==0)
    assert(isHomogeneous D.dd_0)
    assert(degree D=={3})
    assert(prune homology D==cokernel matrix{{x*y,x^2}})
///

TEST /// --test minimizeDM 
    S = QQ[x,y]
    m = matrix{{0,x,y,1},{0,0,0,-y},{0,0,0,x},{0,0,0,0}}
    phi = map(S^{0,1,1,2}, S^{0,1,1,2} ,m, Degree=>2)
    D = differentialModule phi
    M = minimizeDM D
    assert(M.dd_1^2==0)
    assert(isHomogeneous M.dd_0)
    assert(degrees M_0=={{-1},{-1}})
    assert(degree M=={2})
///

TEST /// --test minimizeDM 2
    S = QQ[x,y]
    m = matrix{{0,x^2,x*y,1},{0,0,0,-y},{0,0,0,x},{0,0,0,0}}
    phi = map(S^{0,1,1,3}, S^{0,1,1,3} ,m, Degree=>3)
    D = differentialModule phi
    M = minimizeDM D
    delM = map(S^{1,1},S^{1,1},matrix{{x^2*y,x*y^2},{-x^3,-x^2*y}},Degree=>3)
    assert(differential M==delM)
///

TEST /// --test resDM
    S = QQ[x,y]
    m = matrix{{x*y,y^2},{-x^2,-x*y}}
    phi = map(S^2, S^2, m, Degree=>2)
    D = differentialModule phi
    F = resDM D
    del = map(S^{-1,0,0,1},S^{-1,0,0,1},matrix{{0,x,y,1},{0,0,0,y},{0,0,0,-x},{0,0,0,0}}, Degree=>2)
    assert(F.dd_0^2==0)
    assert(isHomogeneous F.dd_0)
    assert(degree F=={2})
    assert(differential F==del)
///


TEST /// --test resKC
    S = QQ[x,y]
    m = matrix{{x*y,y^2},{-x^2,-x*y}}
    phi = map(S^2, S^2, m, Degree=>2)
    D = differentialModule phi
    F = resKC D
    del = map(S^{-1,0,0,1},S^{-1,0,0,1},matrix{{0,-y,-x,-1},{0,0,0,x},{0,0,0,-y},{0,0,0,0}}, Degree=>2)
    assert(F.dd_0^2==0)
    assert(isHomogeneous F.dd_0)
    assert(degree F=={2})
    assert(differential F==del)
///

TEST /// --test foldComplex
    S = QQ[x,y,z]
    K = koszul vars S
    F0 = foldComplex(K,0)
    F1 = foldComplex(K,1)
    F4 = foldComplex(K,4)
    assert(isHomogeneous differential F1)
    assert(degree F0=={0})
    assert(degree F1=={1})
    assert(degree F4=={4})
    assert(F4.dd_0^2==0)
///

TEST /// --test unfold
    S = QQ[x,y]
    phi = map(S^{1,1},S^{1,1},matrix{{x^2*y,x*y^2},{-x^3,-x^2*y}},Degree=>3)
    D = differentialModule phi
    C = unfold(D,-2,2)
    assert(C.dd_0==D.dd_0)
    assert(C.dd_1==C.dd_0)
    assert(degree C.dd_0=={3})
    assert(C_-2==C_3)
///

TEST /// --resDM 3 vars
    S = QQ[x,y,z]
    phi = map(S^4, S^4, matrix{{x*y,y^2,z,0},{-x^2,-x*y,0,z},{0,0,-x*y,-y^2},{0,0,x^2,x*y}})
    D = differentialModule phi
    F = resDM D    
    assert(F.dd_0^2==0)
///
*-
end;

restart;
load("DifferentialModules.m2");
R = ZZ/101[x,y]
f = map(coker vars R, coker vars R, 0)
D = differentialModule(f)
F = resMinFlag(D, 3)
F.dd


D1 = differentialModule(M)
F = resKC(D1)
F1 = resMinFlag(D1)


S = ZZ/101[x,y]
C = koszul matrix{{x,y}}
D = foldComplex(C, 0)
minFlagOneStep(D)

-- S = ZZ/101[x,y]
-- C = koszul matrix{{x,y}}
-- D = foldComplex(C, 2)
-- minFlagOneStep(D)

-- S = QQ[x,y,z]
-- phi = map(S^4, S^4, matrix{{x*y,y^2,z,0},{-x^2,-x*y,0,z},{0,0,-x*y,-y^2},{0,0,x^2,x*y}})
-- D = differentialModule phi
-- minFlagOneStep(D)
-- F = resDM D  

-- S = QQ[x,y]
-- f = map(S^{0,-1,-1,-2}, S^{0,-1,-1,-2}, matrix{{0,x,y,0},{0,0,0,-y},{0,0,0,x},{0,0,0,0}})

-- g = map(S^{-1,-2,-2,-3}, S^{-1,-2,-2,-3}, matrix{{0,x,y,0},{0,0,0,-y},{0,0,0,x},{0,0,0,0}})

-- h = map(S^{0,-1,-1,-2,-1,-2,-2,-3},S^{0,-1,-1,-2,-1,-2,-2,-3}, f ++ g)

-- D = differentialModule(h)
-- C1 = resMinFlag(D, 4)
-- C2 = resKC(D)
-- C1.dd_1^2
-- trim HH_0(C2) == trim HH_0(D)

-- minFlagOneStep(D)

-- killingCyclesOneStep(D)

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
S = ring weightedProjectiveSpace {1,1,1,2,2};
I = minors(2, matrix{{x_0, x_1, x_2^2, x_3}, {x_1, x_2, x_3, x_4}});
M = Ext^3(module S/I, S^{{-7}});
L = stronglyLinearStrand M
betti L
L.dd
S = ring hirzebruchSurface 3;
M = coker vars S;
stronglyLinearStrand(M)--should give Koszul complex, and it does.
M = coker matrix{{x_1, x_2^2}}
res M
L = stronglyLinearStrand(M)
L.dd
loadPackage "TateOnProducts"
M = coker(N.dd_1)
ring M
S = ring M
A = coker map(S^{-1} ++ S^{-1} ++ S^{-1}, S^{-2} ++ S^{-2} ++ S^{-2} ++ S^{-2} ++ S^{-3} ++ S^{-3}, matrix{{x_0, 0, x_1, 0, x_3, 0}, {0,x_0,-x_2,x_1,-x_4,x_3}, {-x_2,-x_1,0,-x_2,0,-x_4}})
isIsomorphic(M, A)

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
end;


load "MultigradedBGG.m2"


