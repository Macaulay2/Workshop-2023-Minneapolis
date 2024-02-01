-*
newPackage("DifferentialModules",
    Version => "1.1",
    Date => "5 June 2023",
    Headline => "Computing Free Resolutions of Differential Modules",
    Authors => {
	{Name => "Maya Banks",             Email => "mdbanks@wisc.edu",      HomePage => "https://sites.google.com/wisc.edu/mayabanks" },
    {Name => "Michael K. Brown",       Email => "mkb0096@auburn.edu",    HomePage => "http://webhome.auburn.edu/~mkb0096/" },
	{Name => "Daniel Erman",           Email => "erman@wisc.edu",        HomePage => "https://people.math.wisc.edu/~erman/" },
	{Name => "Tara Gomes",             Email => "gomes072@umn.edu",      HomePage => "Fill in" },
	{Name => "Pouya Layeghi",          Email => "layeg001@umn.edu",      HomePage => "Fill in" },
	{Name => "Prashanth Sridhar",      Email => "pzs0094@auburn.ed",     HomePage => "https://sites.google.com/view/prashanthsridhar/home" },
	{Name => "Andrew Tawfeek",         Email => "atawfeek@uw.edu",       HomePage => "https://www.atawfeek.com/" },
	{Name => "Eduardo Torres Davila",  Email => "torre680@umn.edu",      HomePage => "https://etdavila10.github.io/" },
	{Name => "Jay Yang",               Email => "jayy@wustl.edu",        HomePage => "https://www.math.wustl.edu/~jayy/" },
	{Name => "Sasha Zotine",           Email => "18az45@qiueensu.ca",    HomePage => "Fill in" }

	    },
  DebuggingMode => false
  )

export {
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
    assert(#d == 1 and first d == 0);
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
resKC = method();
resKC(DifferentialModule,ZZ) := (D,k)->(
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
resKC(DifferentialModule) := (D)->(
    k := dim ring D + 1;
    d := degree D;
    R := ring D;
    s := numgens D_0;
    scan(k,i-> D = killingCyclesOneStep(D));
    t := numgens D_0;
    newDiff := submatrix(D.dd_1, toList(s..t-1),toList(s..t-1));
    differentialModule (chainComplex(newDiff**R^{-d},newDiff**R^{-d})[1])
    )


--  This code uses the Cartan-Eilenberg construction to build a "free resolution"
--  of an arbitrary differential module.
--  It could use some work.
--  SHOULD WE ADD "prune" when resolviing kernel and image?

resDM = method(TypicalValue => DifferentialModule,
            Options => {
    	    LengthLimit => 3
    	    });

--MAYA: changed to be compatible with new degree convention
resDM(DifferentialModule) := opts -> (D) ->(
    n := opts.LengthLimit;
    d := degree D;
    cyc := res(kernel D,LengthLimit =>n);
    bou := res(image D,LengthLimit =>n);
    R := ring D;
    m := max(length cyc, length bou);
    -- this is code which turns a free complex into a differential module...
    L := apply(length cyc+1,j->(
	    transpose concatMatrices apply(m+1,i->map(cyc_j,cyc_i, if i == j+1 then cyc.dd_i else 0))
	));
    cycDiff := transpose(concatMatrices L);
    cycMod := cyc_0;
    scan(m+1, i-> cycMod = cycMod ++ ((cyc_(i+1))**(R)^{(i+1)*d}));
    L' := apply(length bou+1,j->(
	    transpose concatMatrices apply(m+1,i->map(bou_j,bou_i, if i == j+1 then bou.dd_i else 0))
	));
    bouDiff := transpose(concatMatrices L');
    bouMod := bou_0**(R)^{d};
    scan(m+1, i-> bouMod = bouMod ++ ((bou_(i+1))**(R)^{(i+2)*d}));
    --cycDiff and bouDiff are the differentials of cycle and boundary resolutions, as a diff module
    epsB0 := map(D_0,bou_0, gens image D.dd_0 // gens image D.dd_0);
    epsZ0 := map(D_0,cyc_0,gens kernel D);
    psi := -map(cyc_0,bou_0,gens image D // gens kernel D);
    PSI := {psi};
    scan(m,i-> PSI = PSI|{-PSI_(i)*bou.dd_(i+1) // cyc.dd_(i+1)} );
    tau := map(cyc_0,bou_1, -epsB0*bou.dd_1 // epsZ0);
    dummy := map(bou_0,0,0);
    TAU := {dummy,tau};
    scan(toList(2..m),i-> TAU = TAU|{-TAU_(i-1)*bou.dd_i // cyc.dd_(i-1)});
--I want TAU_i to be the map from bou_i, so we adjoin a dummy variable.
    K := apply(length cyc + 1,j->(
	    apply(length bou+1,i->(
		    map(cyc_j,bou_i, if j == i then PSI_i else 0)
		    ))));
    bouToCycPSI := transpose concatMatrices apply(K,i-> transpose concatMatrices i);
    K' := apply(length cyc + 1,j->(
	    apply(length bou+1,i->(
		    map(cyc_j,bou_i, if j == i - 1 then TAU_i else 0)
		    ))));
    bouToCycTAU := transpose concatMatrices apply(K',i-> transpose concatMatrices i);
    bouToCyc := bouToCycPSI + bouToCycTAU;
    lastMap := map(bouMod,cycMod**R^{-d},0);
    --MAYA CEdiff := map(cycMod++bouMod,(cycMod++bouMod)**(ring D)^{-d}, (cycDiff|bouToCyc) || (lastMap|bouDiff));
    CEdiff := map(cycMod++bouMod,(cycMod++bouMod), (cycDiff|bouToCyc) || (lastMap|bouDiff), Degree=>d); --maya added
    -- MAYA differentialModule(chainComplex(CEdiff**(ring D)^{d},CEdiff)[1])
    differentialModule(chainComplex(CEdiff,CEdiff)[1]) --maya added
	)

--  same as resDM except it outputs a pair (r,eps)
--  where r is the resDM and epsilon is the map r-->D.
--  mostly used to check that resDM produced actual resolutions.
resDMwMap = method(TypicalValue => DifferentialModule,
            Options => {
    	    LengthLimit => 3
    	    });

--MAYA: Changed to be compatible with new degree convention
resDMwMap(DifferentialModule) := opts -> (D) ->(
    n := opts.LengthLimit;
    d := degree D;
    cyc := res(kernel D,LengthLimit =>n);
    bou := res(image D,LengthLimit =>n);
    R := ring D;
    m := max(length cyc, length bou);
    -- this is code which turns a free complex into a differential module...
    L := apply(length cyc+1,j->(
	    transpose concatMatrices apply(m+1,i->map(cyc_j,cyc_i, if i == j+1 then cyc.dd_i else 0))
	));
    cycDiff := transpose(concatMatrices L);
    cycMod := cyc_0;
    scan(m+1, i-> cycMod = cycMod ++ ((cyc_(i+1))**(R)^{(i+1)*d}));
    L' := apply(length bou+1,j->(
	    transpose concatMatrices apply(m+1,i->map(bou_j,bou_i, if i == j+1 then bou.dd_i else 0))
	));
    bouDiff := transpose(concatMatrices L');
    bouMod := bou_0**(R)^{d};
    scan(m+1, i-> bouMod = bouMod ++ ((bou_(i+1))**(R)^{(i+2)*d}));
    --cycDiff and bouDiff are the differentials of cycle and boundary resolutions, as a diff module
    epsB0 := map(D_0,bou_0, gens image D.dd_0 // gens image D.dd_0);
    epsZ0 := map(D_0,cyc_0,gens kernel D);
    psi := -map(cyc_0,bou_0,gens image D // gens kernel D);
    PSI := {psi};
    scan(m,i-> PSI = PSI|{-PSI_(i)*bou.dd_(i+1) // cyc.dd_(i+1)} );
    tau := map(cyc_0,bou_1, -epsB0*bou.dd_1 // epsZ0);
    dummy := map(bou_0,0,0);
    TAU := {dummy,tau};
    scan(toList(2..m),i-> TAU = TAU|{-TAU_(i-1)*bou.dd_i // cyc.dd_(i-1)});
--I want TAU_i to be the map from bou_i, so we adjoin a dummy variable.
    K := apply(length cyc + 1,j->(
	    apply(length bou+1,i->(
		    map(cyc_j,bou_i, if j == i then PSI_i else 0)
		    ))));
    bouToCycPSI := transpose concatMatrices apply(K,i-> transpose concatMatrices i);
    K' := apply(length cyc + 1,j->(
	    apply(length bou+1,i->(
		    map(cyc_j,bou_i, if j == i - 1 then TAU_i else 0)
		    ))));
    bouToCycTAU := transpose concatMatrices apply(K',i-> transpose concatMatrices i);
    bouToCyc := bouToCycPSI + bouToCycTAU;
    lastMap := map(bouMod,cycMod**R^{-d},0);
    CEdiff := map(cycMod++bouMod,(cycMod++bouMod), (cycDiff|bouToCyc) || (lastMap|bouDiff), Degree=>d);-- maya added
    --CEdiff := map(cycMod++bouMod,(cycMod++bouMod)**(ring D)^{-d}, (cycDiff|bouToCyc) || (lastMap|bouDiff));
    --RD := differentialModule(chainComplex(CEdiff**(ring D)^{d},CEdiff)[1]);
    RD := differentialModule(chainComplex(CEdiff,CEdiff)[1]);
    epsBou := matrix{ {epsB0} | apply(length bou, i-> map(D_0,bou_(i+1),0))    };
    epsCyc := matrix{ {epsZ0} | apply(length cyc, i-> map(D_0,cyc_(i+1),0))    };
    eps := epsCyc|epsBou;
    (RD,eps)
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

-- restart;
-- load("DifferentialModules.m2");

--- run examples below here

-- R = ZZ/101[x,y]
-- A = map(R^2, R^2, matrix{{x*y, -x^2}, {y^2, -x*y}}, Degree => 2)
-- D1 = differentialModule(A)
-- F1 = resMinFlag(D1, 4)

-- F = resKC(D1)


-- S = ZZ/101[x,y]
-- C = koszul matrix{{x,y}}
-- D = foldComplex(C, 0)
-- minFlagOneStep(D)

-- S = ZZ/101[x,y]
-- C = koszul matrix{{x,y}}
-- D = foldComplex(C, 2)
-- minFlagOneStep(D)

-- S = QQ[x,y,z]
-- phi = map(S^4, S^4, matrix{{x*y,y^2,z,0},{-x^2,-x*y,0,z},{0,0,-x*y,-y^2},{0,0,x^2,x*y}})
-- D = differentialModule phi
-- << degree(D) << endl;
-- minFlagOneStep(D)
-- F = resMinFlag

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
