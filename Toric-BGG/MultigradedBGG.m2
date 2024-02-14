newPackage("MultigradedBGG",
    Version => "1.1",
    Date => "5 June 2023",
    Headline => "the multigraded BGG correspondence and differential modules",
    Authors => {
	{Name => "Maya Banks",         	     Email => "mdbanks@wisc.edu",      HomePage => "https://sites.google.com/wisc.edu/mayabanks" },
        {Name => "Michael K. Brown",         Email => "mkb0096@auburn.edu",    HomePage => "http://webhome.auburn.edu/~mkb0096/" },
	{Name => "Daniel Erman",    	     Email => "erman@wisc.edu",        HomePage => "https://people.math.wisc.edu/~erman/" },
	{Name => "Tara Gomes",	    	     Email => "gomes072@umn.edu",      HomePage => "Fill in" },
	{Name => "Prashanth Sridhar",	     Email => "pzs0094@auburn.edu",    HomePage => "https://sites.google.com/view/prashanthsridhar/home"},
	{Name => "Eduardo Torres Davila",    Email => "torre680@umn.edu",      HomePage => "https://etdavila10.github.io/" },
	{Name => "Sasha	Zotine",    	     Email => "18az45@queensu.ca",     HomePage => "https://sites.google.com/view/szotine/home" }
    },
    PackageExports => {"NormalToricVarieties", "Complexes"},
    DebuggingMode => true
  )

export {
    --methods
    "dualRingToric",
    "toricRR",
    "toricLL",
    "stronglyLinearStrand",
    "unfold",
    "DifferentialModule",
    "differentialModule",
    "foldComplex",
    "resDM",    
    "minimizeDM",
    "differential",
    --symbols
    "SkewVariable"
    }

--------------------------------------------------
--- Differential Modules
--------------------------------------------------
DifferentialModule = new Type of Complex
DifferentialModule.synonym = "differential module"

differentialModule = method(TypicalValue => DifferentialModule)
differentialModule Complex := C -> (
    assert(C == naiveTruncation(C, -1, 1)) --checks if C is concentrated in homological degrees -1, 0, 1
    --add check that the two differentials
    --are the same matrix (say of degree d, if there 
    --is a grading), this matrix squares to 0, C_1 = C_0(-d), and C_{-1} = C_0(d).
  new DifferentialModule from C);
differentialModule Matrix := phi -> (
    --check if the source and target are the same up to a twist
    if phi^2 != 0 then error "The differential does not square to zero.";
    R := ring phi;
    if target phi != source phi then error "source and target of map are not the same"; 
    new DifferentialModule from (complex({-phi,-phi})[1]));--the maps are negated to cancel out the sign introduced by the shift.

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
degree DifferentialModule := ZZ => (D -> degree D.dd_1); 
differential = method();
differential DifferentialModule := Matrix=> (D->D.dd_1);
kernel DifferentialModule := Module => opts -> (D -> kernel D.dd_0); 
image DifferentialModule := Module => (D -> image D.dd_1); 
homology DifferentialModule := Module => opts -> (D -> HH_0 D);

isFreeModule(DifferentialModule) := D ->(
    isFreeModule module D
    )

unfold = method();
--Input:  a differential module and a pair of integers low and high
--Output:  the unfolded chain complex of the differential module, in homological degrees
--         low through high.
unfold(DifferentialModule,ZZ,ZZ) := Complex => (D,low,high)->(
    L := toList(low..high);
    d := degree D;
    R := ring D;
    phi := differential D;
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
    assert(first d == 0);
    R := ring D;
    s := numgens D_0;
    scan(k,i-> D = minFlagOneStep(D));
    t := numgens D_0;
    newDiff := submatrix(D.dd_1, toList(s..t-1),toList(s..t-1));
    differentialModule (complex({-newDiff**R^{-d},-newDiff**R^{-d}})[1])
)

        
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
    newDiff := matrix{{D.dd_1,psi},{map(G_0**R^{d},D_1,0, Degree=>d),map(G_0**R^{d},G_0**R^{d},0, Degree=>d)}}; 
    assert (newDiff*(newDiff) == 0);
    differentialModule newDiff
    )
--killing cycles resolution
--Input:  a free(?) differential module D and an integer k
--Output:  the killing cycles resolution of D after k steps.
-- NOTE:  until recently this produced the CONE over the resolution.
--        but I'm interested in the resolution itself so the last few lines change that
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
    differentialModule (complex({-newDiff**R^{-d},-newDiff**R^{-d}})[1])
    )

--same as above, but default value of k 
resDM(DifferentialModule) := (D)->(
    k := dim ring D + 1;
    d := degree D;
    R := ring D;
    s := numgens D_0;
    scan(k,i-> D = killingCyclesOneStep(D));
    t := numgens D_0;
    newDiff := submatrix(D.dd_1, toList(s..t-1),toList(s..t-1));
    differentialModule (complex({-newDiff**R^{-d},-newDiff**R^{-d}})[1])
    )

--  Subroutines and routines to produce the minimal part of a matrix.
--  NEEDS TO BE A SQUARE MATRIX
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
minimizeDM = method();
minimizeDM(DifferentialModule) := r ->(
    R := ring r;
    d := degree r;
    A := minimizeDiff(r.dd_1);
    degA := map(target A, source A, A, Degree=>d);
    differentialModule (complex({-degA,-degA})[1])
    )

---
---

--  Input:  a free complex F and a degree d
--  Output: the corresponding free differential module of degree da
foldComplex = method();
foldComplex(Complex,ZZ) := DifferentialModule => (F,d)->(
    R := ring F;
    L := apply(length F+1,j->(
	    --sasha: i removed the concatMatrices method since 'matrix {_}' just does the same thing
	    transpose matrix {apply(length F+1,i->map(F_j,F_i, if i == j+1 then F.dd_i else 0))}
	));
    FDiff := transpose(matrix {L});
    FMod := F_0;
    scan(length F+1, i-> FMod = FMod ++ ((F_(i+1))**(R)^{(i+1)*d}));
    degFDiff := map(FMod,FMod,FDiff, Degree=>d); 
    differentialModule(complex({-degFDiff,-degFDiff})[1]) 
    )

--------------------------------------------------
--- Toric BGG
--------------------------------------------------

-- Exported method for dualizing an algebra.
-- Input:   a polynomial ring OR an exterior algebra
-- Output:  the Koszul dual exterior algebra, but with an additional 
--          ZZ-degree, the ``standard grading'' where elements of \bigwedge^k
--          have degree k.
dualRingToric = method(Options => {
    	    Variable        => getSymbol "x",
	    SkewVariable => getSymbol "e"});
dualRingToric PolynomialRing := opts -> S ->(
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

-- Exported method for computing the multigraded BGG functor of a module over
-- a polynomial ring.
-- Input:   M a (multi)-graded S-module.
--          L a list of degrees.
-- Output:  The differential module RR(M) in degrees from L,
--          presented as a complex in homological degrees -1, 0 ,1
--          and with the same differential in both spots.
toricRR = method();
toricRR(Module,List) := (N,L) ->(
    M := coker presentation N;
    S := ring M;
    if not isCommutative S then error "--base ring is not commutative";
    E := dualRingToric S;
    -- pick out the generators in the selected degree range.
    f0 := matrix {for d in unique L list gens image basis(d,M)};
    -- this is the degree of the anti-canonical bundle, which
    -- we need to twist by when we take the koszul dual.
    wEtwist := append(-sum degrees S, -numgens S);
    -- here are the degrees of the generators we will have in the dual.
    f0degs := apply(degrees source f0, d -> (-d | {0}) + wEtwist);
    if #f0degs == 0 then complex map(E^0, E^0, 0) else (	
    	relationsM := presentation M;
    	SE := S**E;
    	tr := sum(dim S, i-> SE_(dim S+i)*SE_i);
    	newf0 := sub(f0,SE)*tr;
    	relationsMinSE := sub(relationsM,SE);
    	newf0 = newf0 % relationsMinSE;
    	newg := contract(transpose sub(f0,SE), newf0);
	-- null input is to make the map have the correct degrees to be homogeneous.
	newg = transpose map(source newg, , transpose newg);
    	g' := sub(newg,E);
	del := map(E^f0degs,E^f0degs, -g', Degree => degree 1_S | {-1});
	differentialModule(complex({del, del})[1])
    	)
    )

-- If there is no input of L, we make a simple choice based on some information on M.
toricRR Module := M -> (
    L := if length M < infinity then unique flatten degrees basis M
    else join(degrees M, apply(degrees ring M, d -> (degrees M)_0 + d));
    toricRR(M,L)
    )

-- Exported method for computing the multigraded BGG functor of a module over the exterior algebra.
-- Input: N a (multi)-graded E-module.
-- Output: a complex of modules over the polynomial ring.
-- Caveat: Assumes N is finitely generated.
-- Caveat 2:  arrows of toricLL(N) correspond to exterior multiplication (not contraction)
toricLL = method();
toricLL Module := N -> (
    E := ring N;
    if not isSkewCommutative E then error "ring N is not skew commutative";
    S := dualRingToric E;
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
    if #homDegs == 1 then (complex {map(S^0,FF#0,0)})[1] else (
	complex apply(drop(homDegs,-1), i-> map(FF#i,FF#(i+1), (-1)^(homDegs#0)*differential_(inds#(i+1))^(inds#(i))))[-homDegs#0]
	)
    )

-- Exported method for computing a strongly linear strand for a module over a polynomial ring.
-- Input: M a (multi)-graded module over a polynomial ring
-- Output: a Complex, the strongly linear strand of the minimal
--     	   free resolution of M (in the sense of the paper "Linear strands  
--    	   of multigraded free resolutions" by Brown-Erman).
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


beginDocumentation()

--------------------------------------------------
--- Differential Modules
--------------------------------------------------

doc ///
   Key 
      MultigradedBGG
   Headline 
      Package for working with Multigraded BGG and Differential Modules
   Description
    Text
      todo    
///


doc ///
   Key 
    differentialModule
    (differentialModule, Matrix)
   Headline
    converts a square zero matrix into a differential module
   Usage
    differentialModule(f)
   Inputs
    f : Matrix
        representing a module map with the same source and target
   Outputs
    : DifferentialModule 
   Description
    Text
      Given a module $f: M\to M$ of degree a this creates a degree a differential module from
      f represented as a 3-term chain complex in degree -1, 0, 1. If you want a nonzero, 
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
    : Complex
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
    (foldComplex,Complex,ZZ)
   Headline
    converts a chain complex into a differential module
   Usage
    foldComplex(C)
   Inputs
    C: Complex
    d: ZZ
   Outputs
    : DifferentialModule
   Description
    Text
      Given a chain complex C and integer d it creates the corresponding
      (flag) differential module of degree d.
    Example
      R = QQ[x,y];
      C = complex res ideal(x,y)
      D = foldComplex(C,0);
      D.dd_1
///

--todo: since resKC has been combined with resDM, which of these doc nodes is needed? can we cut
-- some of the above into this node?
doc ///
   Key 
    resDM
    (resDM,DifferentialModule)
    (resDM,DifferentialModule,ZZ)
   Headline
    uses a "killing cycles"-style construction to find a free resolution of a differential module
   Usage
    resDM(D)
    resDM(D,k)
   Inputs
    D: DifferentialModule
    k: ZZ
   Outputs
    : DifferentialModule
   Description
    Text
      Given a differential module D it creates a free flag resolution of D, using a killing cycles
      construction.  Because of issues with adding options, there are two choices.  The default
      resDM(D) runs the algorithm for the number of steps determined by the dimension of the ambient ring.
      resDM(D,k) for k steps.
    Example
      R = QQ[x,y];
      M = R^1/ideal(x^2,y^2);
      phi = map(M,M,x*y, Degree=>2);
      D = differentialModule phi;
      r = resDM(D)
      r.dd_1
    Text
      The default algorithm runs for dim R + 1 steps, again because I couldn't figure out the options.
      Adding the number of steps as a second argument is like adding a LengthLimit.
    Example
      R = QQ[x]/(x^3);
      phi = map(R^1,R^1,x^2, Degree=>2);
      D = differentialModule phi;
      r = resDM(D)
      r.dd_1      
      r = resDM(D,6)
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
      r = resDM(D)
      r.dd_1
      mr = minimizeDM(r)
      mr.dd_1   
///

--------------------------------------------------
--- Toric BGG
--------------------------------------------------

doc ///
   Key 
    dualRingToric
    (dualRingToric, PolynomialRing)
    [(dualRingToric, PolynomialRing),Variable, SkewVariable]
   Headline
    computes the Koszul dual of a polynomial ring/exterior algebra
   Usage
    dualRingToric R
   Inputs
    R : PolynomialRing
        either a standard polynomial ring, or an exterior algebra
    Variable : Symbol
    SkewVariable : Symbol
   Outputs
    : PolynomialRing
      which is the Koszul dual of R
   Description
    Text
      This method computes the Koszul dual of a polynomial ring or exterior algebra. In
      particular, if $S = k[x_0, \ldots, x_n]$ is a $\mathbb{Z}^m$-graded ring for some $m$, 
      and $\operatorname{deg}(x_i) = d_i$, then the output of dualRingToric is the 
      $\mathbb{Z}^{m+1}$-graded exterior algebra on variables $e_0, \ldots, e_n$ with degrees 
      $(-d_i, -1)$.
    Example
      R = ring(hirzebruchSurface(2, Variable => y))
      E = dualRingToric(R, SkewVariable => f)
    Text
      On the other hand, if $E$ is a $\mathbb{Z}^{m+1}$-graded exterior algebra on $n+1$ 
      variables $e_0, \ldots, e_n$ with $\operatorname{deg}(e_i) = (-d_i, -1)$, then 
      dualRingToric E is the $\mathbb{Z}^m$-graded polynomial ring $k[x_0, \ldots, x_n]$ with 
      $\operatorname{deg}(x_i) = d_i$.
    Example
      RY = dualRingToric E
      degrees RY == degrees R
    Text
      This method preserves the coefficient ring of the input ring.
    Example
      S = ZZ/101[x,y,z, Degrees => {1,1,2}]
      E' = dualRingToric S
      coefficientRing E' === coefficientRing S
///

doc ///
   Key 
    toricRR
    (toricRR, Module)
    (toricRR, Module, List)
   Headline
    computes the BGG functor of a module over a multigraded polynomial ring. 
   Usage
    toricRR M
    toricRR(M,L)
   Inputs
    M : Module
        a module over a multigraded polynomial ring
    L : List 
        a list of multidegrees of the polynomial ring
   Outputs
    : DifferentialModule
      a quotient of the differential module obtained by applying the multigraded BGG functor R 
      to M. The size of this quotient is determined by the list L. 
   Description
    Text
       todo
    Example
       a
///

doc ///
   Key 
    toricLL
    (toricLL, Module)
   Headline
    computes the BGG functor of a module over a multigraded exterior algebra
   Usage
    toricLL N
   Inputs
    N : Module
        a module over a multigraded exterior algebra
   Outputs
    : Complex
      a complex of modules over the polynomial ring which is the output of the BGG functor
   Description
    Text
      todo 
    Example
      R = ring(hirzebruchSurface(2, Variable => y))
///

doc ///
   Key 
    stronglyLinearStrand
    (stronglyLinearStrand, Module)
   Headline
    computes the strongly linear strand of the minimal free resolution of a finitely generated 
    multigraded module over a polynomial ring, provided the module is generated in a single degree
   Usage
    stronglyLinearStrand N
   Inputs
    M : Module
        a finitely generated module over a multigraded polynomial ring that is generated in a
	single degree
   Outputs
    : Complex
      the strongly linear strand of the minimal free resolution of M
   Description
    Text
      The strongly linear strand is defined in the paper "Linear strands of multigraded free 
      resolutions" by Brown-Erman. It is, roughly speaking,the largest subcomplex of the minimal 
      free resolution of M that is linear, in the sense that its differentials are matrices of 
      linear forms. The method we use for computing the strongly linear strand uses BGG, mirroring
      Corollary 7.11 in Eisenbud's textbook "The geometry of syzygies".
    Example
      R = ring(hirzebruchSurface(2, Variable => y))
///




--- TESTS

--------------------------------------------------
--- Differential Modules
--------------------------------------------------

--Michael: I changed some signs in the differential modules methods that may affect these tests. Some signs
--may need to be flipped. 

-- Constructor tests
TEST /// 
    S = QQ[x,y]
    m = matrix{{0,x,y,1},{0,0,0,-y},{0,0,0,x},{0,0,0,0}}
    phi = map(S^{0,1,1,2}, S^{0,1,1,2} ,m, Degree=>2)
    D = differentialModule phi
    assert(D.dd_0^2==0)
    assert(isHomogeneous D.dd_0)
    assert(degree D=={2})
    assert(prune homology D==cokernel matrix{{x,y}})
///

TEST ///
    S = QQ[x,y]
    m = matrix{{0,x^2,x*y,1},{0,0,0,-y},{0,0,0,x},{0,0,0,0}}
    phi = map(S^{0,1,1,3}, S^{0,1,1,3} ,m, Degree=>3)
    D = differentialModule phi
    assert(D.dd_0^2==0)
    assert(isHomogeneous D.dd_0)
    assert(degree D=={3})
    assert(prune homology D==cokernel matrix{{x*y,x^2}})
///

-- Testing minimizeDM
TEST /// 
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

TEST ///
    S = QQ[x,y]
    m = matrix{{0,x^2,x*y,1},{0,0,0,-y},{0,0,0,x},{0,0,0,0}}
    phi = map(S^{0,1,1,3}, S^{0,1,1,3} ,m, Degree=>3)
    D = differentialModule phi
    M = minimizeDM D
    delM = map(S^{1,1},S^{1,1},matrix{{x^2*y,x*y^2},{-x^3,-x^2*y}},Degree=>3)
    assert(differential M==delM)
///

-- Testing resDM
TEST ///
    S = QQ[x,y]
    m = matrix{{x*y,y^2},{-x^2,-x*y}}
    phi = map(S^2, S^2, m, Degree=>2)
    D = differentialModule phi
    F = resDM D
    del = map(S^{-1,0,0,1},S^{-1,0,0,1},matrix{{0,-y,-x,-1},{0,0,0,x},{0,0,0,-y},{0,0,0,0}}, Degree=>2)
    assert(F.dd_0^2==0)
    assert(isHomogeneous F.dd_0)
    assert(degree F=={2})
    assert(differential F==del)
///

TEST ///
    S = QQ[x,y,z]
    phi = map(S^4, S^4, matrix{{x*y,y^2,z,0},{-x^2,-x*y,0,z},{0,0,-x*y,-y^2},{0,0,x^2,x*y}})
    D = differentialModule phi
    F = resDM D    
    assert(F.dd_0^2==0)
///

-- Testing foldComplex
TEST ///
    S = QQ[x,y,z]
    K = koszulComplex vars S
    F0 = foldComplex(K,0)
    F1 = foldComplex(K,1)
    F4 = foldComplex(K,4)
    assert(isHomogeneous differential F1)
    assert(degree F0=={0})
    assert(degree F1=={1})
    assert(degree F4=={4})
    assert(F4.dd_0^2==0)
///

-- Testing unfold
TEST ///
    S = QQ[x,y]
    phi = map(S^{1,1},S^{1,1},matrix{{x^2*y,x*y^2},{-x^3,-x^2*y}},Degree=>3)
    D = differentialModule phi
    C = unfold(D,-2,2)
    assert(C.dd_0==D.dd_0)
    assert(C.dd_1==C.dd_0)
    assert(degree C.dd_0=={3})
    assert(C_-2==C_3)
///

--------------------------------------------------
--- Toric BGG
--------------------------------------------------

-- Testing dualRingToric
TEST ///
S = ring(hirzebruchSurface(2, Variable => y));
E = dualRingToric(S,SkewVariable => f);
SY = dualRingToric(E);
assert(degrees SY == degrees S)
///

-- Testing toricRR
TEST ///
S = ring hirzebruchSurface 3;
M = coker matrix{{x_0}};
L = {{0,0}, {1,0}};
D = toricRR(M, L);
assert(degree D == {0,0,-1})
E = ring D;
f = map(E^{{1, -2, -4}} ++ E^{{0, -2, -4}}, E^{{1, -2, -4}} ++ E^{{0, -2, -4}}, matrix{{0, 0}, {e_2, 0}});
assert(D.dd_0 == f)
L = {{0,0}, {1,0}, {-3, 1}, {0,1}, {2,0}};
D = toricRR(M,L);
assert(degree D == {0,0,-1})
assert(D.dd^2 == 0)
--this takes several seconds, so we won't include it
--M = coker random(S^2, S^{3:{-3,-2}});
--D = toricRR M
--assert(degree D == {0,0,-1})
--assert(D.dd^2 == 0)
///

TEST ///
S = ring weightedProjectiveSpace {1,1,2}
M = coker random(S^2, S^{3:{-5}});
D = toricRR M
assert(degree D == {0,-1})
assert(D.dd^2 == 0)
///

--Testing toricLL
TEST ///
S = ring hirzebruchSurface 3;
E = dualRingToric S;
C = toricLL(E^1)
assert (isHomogeneous C)
assert ((C.dd)^2 == 0)
--The following 5 assertions check that C is isomorphic to the 
--Koszul complex on x_0, ..., x_3 (up to a twist and shift)
assert (HH_(0) C == 0)
assert (HH_(-1) C == 0)
assert (HH_(-2) C == 0)
assert (HH_(-3) C == 0)
assert (minors(1, C.dd_(-3)) == ideal vars ring C)
N = coker vars E
C = toricLL(N)
assert(rank C_0 == 1)
assert(C_-1 == 0)
assert(C_1 == 0)
///

TEST ///
S = ring hirzebruchSurface 3;
E = dualRingToric S;
-- cyclic but non-free module:
C = toricLL(coker matrix{{e_0, e_1}})
assert(C == koszulComplex {x_2, x_3})
///

TEST ///
E = dualRingToric (ZZ/101[x_0, x_1, Degrees => {1, 2}]);
--non-cyclic module
N = module ideal(e_0, e_1);
C = toricLL(N);
S = ring C;
f = map(S^{{3}}, S^{{1}} ++ S^{{2}}, matrix{{x_1, -x_0}});
C' = complex({f})[-2];
assert(C == C')
///

--Testing stronglyLinearStrand 
TEST ///
S = ring weightedProjectiveSpace {1,1,1,2,2};
I = minors(2, matrix{{x_0, x_1, x_2^2, x_3}, {x_1, x_2, x_3, x_4}});
M = Ext^3(module S/I, S^{{-7}});
L = stronglyLinearStrand M;
loadPackage "TateOnProducts"
M = coker(L.dd_1)
S = ring M
A = coker map(S^{-1} ++ S^{-1} ++ S^{-1}, S^{-2} ++ S^{-2} ++ S^{-2} ++ S^{-2} ++ S^{-3} ++ S^{-3}, matrix{{x_0, 0, x_1, 0, x_3, 0}, {0,x_0,-x_2,x_1,-x_4,x_3}, {-x_2,-x_1,0,-x_2,0,-x_4}})
assert(isIsomorphic(M, A))
assert(HH_2 L == 0)
assert(HH_1 L == 0)
--Aside: M is the canonical module of the coordinate ring of a copy of P^1 embedded in
--the weighted projective space P(1,1,1,2,2). 
///

TEST ///
S = ring hirzebruchSurface 3;
M = coker vars S;
L = stronglyLinearStrand(M)--should give Koszul complex, and it does.
assert(numcols basis HH_0 L == 1)
assert(HH_1 L == 0)
assert(HH_2 L == 0)
assert(HH_3 L == 0)
assert(HH_4 L == 0)
M = coker matrix{{x_1, x_2^2}}
L = stronglyLinearStrand(M)
assert(L == koszulComplex {x_1})
///

