needsPackage "TestIdeals"

needsPackage "Saturation"

----------------------------------------------------------------------------------------------
-- Auxiliary Functions
----------------------------------------------------------------------------------------------

FrobeniusFunctor = new Type of MethodFunction

FModule = new Type of HashTable

ModuleClass = new Type of Module

GeneratingMorphism = new Type of Morphism

----------------------------------------------------------------------------------------------
-- ModuleClass 
----------------------------------------------------------------------------------------------

-- Module classes are equal if they have the same minimal presentation
ModuleClass == ModuleClass := ( M, N ) -> 
( 
    M1 := new Module from minimalPresentation M; 
    -- because sometimes minimalPresentation returns ModuleClass, causing infinte recursion
    N1 := new Module from minimalPresentation N;
    M1 == N1
)

moduleClass = method()

moduleClass Module := ModuleClass => M -> new ModuleClass from M

----------------------------------------------------------------------------------------------
-- Frobenius Functor
----------------------------------------------------------------------------------------------

FF = method()

-- TODO: need to check if ring is a polynomial ring (or regular) 

FF ( ZZ, Module ) := Module => ( e, M ) -> 
(
    if isFreeModule M or M == 0 then return M; 
    R := ring M;
    p := char R;
    local N; local degsource; local degtarget;
    if isSubmodule M then
    (   --submodule of a free module
        N = gens M;
        -- adjustments of degrees are needed to keep things homogeneous and maps composable
        degsource = - p^e * ( degrees source N );
        degtarget = - p^e * ( degrees target N );
        image map( R^degtarget, R^degsource, entries frobenius^e N )
    ) 
    else if isQuotientModule M then
    (
        N = relations M;
        -- adjustments of degrees are needed to keep things homogeneous and maps composable
        degsource = - p^e * ( degrees source N );
        degtarget = - p^e * ( degrees target N );
        coker map( R^degtarget, R^degsource, entries frobenius^e N )
    )
    else
    (   -- assuming that this means subquotient
        N = image gens M;
        K := image relations M;
        FF( e, N )/FF( e, K )
    )    
)

FF ( ZZ, Matrix ) := Matrix => ( e, f ) -> 
    map( FF(e, target f), FF(e, source f), entries frobenius^e f )

FF Thing := M -> FF( 1, M )

FF = new FrobeniusFunctor from FF

FrobeniusFunctor ^ ZZ := ( f, n ) -> ( x -> f( n, x ) )

----------------------------------------------------------------------------------------------
-- F-Modules
----------------------------------------------------------------------------------------------

generatingMorphism = method()

-- Given a standard map f, generatingMorphism(f) verifies whether f well defined, and whether
-- f maps a module M to F(M). If so, it returns a GeneratingMorphism identical to f.
generatingMorphism Matrix := Matrix => f ->
(
    if f == 0 then return new GeneratingMorphism from f;
    if not isWellDefined f then 
        error "generatingMorphism: map is not well defined";
    if target f != FF( source f ) then 
        error "generatingMorphism: map does not map a module M to F(M)";
    new GeneratingMorphism from stdMap f
)

generatingMorphism Matrix := GeneratingMorphism => f -> generatingMorphism morphism f

-- Standardizes a map and tries to make a GeneratingMorphism from it.
--generatingMorphism Matrix := GeneratingMorphism => f -> generatingMorphism morphism f
    
makeFModule = method()

makeFModule GeneratingMorphism := FModule => g -> 
    new FModule from { generatingMorphism => g, cache => new CacheTable }

makeFModule Matrix := FModule => g -> makeFModule generatingMorphism  g

--- Compute a generating morphism for H_I^i(R)

localCohomologyExt := ( i, I, R ) -> 
(
    -- TODO: check that I is ideal of R, positive characteristic, etc.
    M := R^1/I;
    f := inducedMap( M, FF M );
    E := Ext^i( f, R^1 );
    makeFModule generatingMorphism E
)

--- Compute a generating morphism for H_I^i(R)
localCohomologyFilter := ( i, I, R ) -> 
(
    filterSeq := randomFilterRegSeq( i, I, R );
    J := ideal filterSeq;
    p := char R;
    u := ( product filterSeq )^( p-1 );
    time K := ascendingIdealEquality( 
        e ->  frobeniusPower( p^e, J ) : ideal( u^(lift( (p^e-1)/(p-1), ZZ )) ) 
    );
    time M1 := saturate( R^1/(frobenius K), I );
    time N1 := saturate( R^1/K, I );
    rtMorphism := generatingMorphism map( M1, N1, u );
    M := makeFModule rtMorphism;
    M#cache#(symbol root) = rtMorphism;
    M
)

localCohomology = method( Options => { Strategy => Ext } )

localCohomology ( ZZ, Ideal, Ring ) := FModule => o -> ( i, I, R ) -> 
(
    -- TODO: check that I is ideal of R, positive characteristic, etc.
    if o.Strategy === Ext then localCohomologyExt( i, I, R )
    else localCohomologyFilter( i, I, R )
)
  
root = method()

root FModule := Module => ( cacheValue symbol root )( M ->
(
    g := M#generatingMorphism;
    K := ker g;
    if K == 0  then 
    (    
        if debugLevel > 1 then 
            print "generatingRoot: generating morphism is already injective"; 
        return source g
    );
    g1 := (FF g)*g;
    K1 := ker g1;
    counter := 1;
    if debugLevel > 1 then 
        print( "generatingRoot: computed kernel #"  | toString counter );
    while K1 != K do
    (
        K = K1;
        g1 = (FF g1)*g;
        K1 = ker g1;
        counter = counter + 1;
        if debugLevel > 1 then 
            print( "generatingRoot: computed kernel #" | toString counter )
    );   
    (source g)/K
))
           
ZZ == FModule := ( n, M ) -> M == n
           
FModule == ZZ := ( M, n ) -> 
(
    if n =!= 0 then error "Attempted to compare an FModule to nonzero integer";
    -- check while generating root == 0
    root( M )  == 0
)

cohomDim = method( Options => { Strategy => Ext } )

cohomDim Ideal := ZZ => o -> I ->
(
    R := ring I;
    n := #(trim I)_*;
    while localCohomology( n, I, R, o) == 0 do n = n-1;
    n
)

associatedPrimes FModule := List => o -> M -> associatedPrimes( root M, o )

----------------------------------------------------------------------------------------------
-- Generating random generating morphisms and FModules
----------------------------------------------------------------------------------------------

randomGeneratingMorphism = method( Options=> { Degree => 1 } )

-- produces a random generating morphism M --> F(M)
randomGeneratingMorphism ModuleClass := generatingMorphism => o -> M ->
(
    H := Hom( M, FF M );
    n := numgens H;
    d := o.Degree;
    H = toList apply( n, i -> homomorphism H_{i} );
    R := ring M;
    coeff := flatten entries random( R^{n:d}, R^1 );    
    generatingMorphism sum( coeff, H, ( x, y ) -> x*y )    
)

----------------------------------------------------------------------------------------------

--- I-filter regular sequences ---
isFilterRegElement = method()

isFilterRegElement (RingElement,Ideal,Module,Module) := Boolean => (x,I,M,N) ->
(
    T:=ker (x*id_(M/N));
    J:=radical(ann(T));
    isSubset(I,J)
)

isFilterRegSeq = method()

isFilterRegSeq (BasicList,Ideal,Module) := Boolean => (L,I,M) ->
(
    if not isSubset(ideal(L),I) then error "isFilterRegSeq: The sequence is not contained in the ideal.";
    l:=#L;
    R:=ring I;
    isFRE = isFilterRegElement(L#0,I,M,ideal(0_R)*M);
    N:=ideal(L#0)*M;
    i:=1;
    while (isFRE and i<l) do (
	isFRE = isFilterRegElement(L#i,I,M,N);
	N=N+ideal(L#i)*M;
	i=i+1;
	);
    isFRE
)

isFilterRegSeq (BasicList,Ideal,Ring) := Boolean => (L,I,R) ->
(
    isFilterRegSeq(L,I,module R)
)

isFilterRegSeq (BasicList,Ideal,Ideal) := Boolean => (L,I,R) ->
(
    isFilterRegSeq(L,I,module R)
)

filterRegSeq = method()

filterRegSeq (ZZ,Ideal,Module) := List => (n,I,M) ->
(
    --ap:= associatedPrimes M;
    --ap:= select(ap,p -> (not isSubset(I,p)))
    
    --Strategy of ordering these generators obtained from Eisenbud-Huneke-Vasconcelos.m2 in the
    ---PrimaryDecomposition package.
    G:= sort(flatten entries mingens I, f -> (sum degree f, #terms f));
    k := coefficientRing ring I;
    f := 1_k;
    
    i:=0;
    while not isFilterRegElement(G#i,I,M,ideal(0_R)*M) do (i=i+1;);
    L:={G#i};
    i=0;
        
    while (#L < n and i<#G) do (
--	 if any(ap, p -> (G#i % p==0))
    	if isFilterRegSeq(append(L,G#i),I,M) then (L=append(L,G#i);
	    print L;
	    );
	i=i+1;
	);
    i=0;
    while (#L<n and i<#G) do (
	for j to (#G-1) do ( 
	    if isFilterRegSeq(append(L,(G#i+f*G#j)),I,M) then (L=append(L,(G#i+f*G#j));
		print L;
		print "hi";
		);
	    if #L==n then return L;
	    );
	i=i+1;
	);
    L
)

----------------------------------------------------------------------------------------------

-- generates a random element of degree deg of the ideal I
randomElementInIdeal = method( Options => { Homogeneous => false } )

randomElementInIdeal ( ZZ, RR, Ideal ) := RingElement => o -> ( deg, density, I ) ->
(   
    if o.Homogeneous or deg == 1 then random( deg, I, Density => density )
    else sum random( toList( 1..deg ), I, Density => density )
)

randomFilterRegSeq = method( 
    Options => 
    { 
        Tries => infinity, 
        Homogeneous => false, 
        MaxDegree => infinity 
    } 
)

randomFilterRegSeq ( ZZ, Ideal, Module ) := List => o -> ( n, I, M ) -> 
(
    L := {};
    G := (trim I)_*;
    minDeg := min apply( G, k -> first degree k );
    minDensity := 0.1;
    R = ring I;
    J := ideal( 0_R ); 
    counter := 0;
    local candidate; local deg; local density;
    while counter < o.Tries and deg < o.MaxDegree and #L < n do
    (
        deg = minDeg + ( counter // 100 );
        density = minDensity + 0.009*( counter % 100); 
        if counter < #G then candidate = G_counter
        else candidate = randomElementInIdeal( deg, density, I, Homogeneous => o.Homogeneous);
        if isFilterRegElement( candidate, I, M, J*M ) then 
        (
            L = append( L, candidate );
            J = J + ideal( candidate )
        );
        counter = counter + 1;
    );
    if #L < n then error "randomFilterRegSeg: could not find a sequence of the desired length; try increasing Tries or MaxDegree";
    L
)

randomFilterRegSeq ( ZZ , Ideal, Ring ) := List => o -> ( n, I, R ) ->
    randomFilterRegSeq( n, I, R^1 )
    


----------------------------------------------------------------------------------------------

-- generates limit closure and lower limit ideal


ascendingIdealEquality = method( Options => { MaxTries => infinity } )
-----this function should eventually check to make sure the chain of ideals is ascending

-----this function specifically finds the FIRST spot that we retain equlaity and returns
-----the union of ideals up to that spot (equivalently the ideal at that spot)
ascendingIdealEquality Function := Ideal => o -> f ->
(
    i := 1;
    isEqual = false;
    while not isEqual and i < o.MaxTries do 
    (
	isEqual = ( f(i) == f(i+1) );
	i = i+1;
    );
    if not isEqual then error "ascendingIdealEquality: Reached maximum limit of tries.";
    f(i)
)

limitClosure = method()
limitClosure BasicList := Ideal => L -> 
(
    if #L == 0 then error "limitClosure: limit closure should be the 0 ideal, but cannot check the ring it lives in since the list is empty.";
    p := char ring L#0;
    ascendingIdealEquality(j -> (frobenius^j(ideal(L)):ideal(product(L,x->x^(p^j-1)))))
)

lowerLimit = method()
lowerLimit ( BasicList, Ring ) := Ideal => ( L, R ) ->
(
    --By convention, I am pretty certain that we return the 0 ideal in the case of an empty
    --list. See Example 6.5 from "Lyubeznik numbers, F-modules, and modules of generalized fractions
    if #L == 0 then return ideal 0_R;
    local LC;
    if #L == 1 then (
	if L#0 == 0_R then return ideal 0_R;
	LC = ideal 0_R;
    )
    else LC = limitClosure drop( L, -1 );
    ascendingIdealEquality( j -> LC : L#(-1)^j )
)

lowerLimit RingElement := Ideal => f -> lowerLimit { f }
    
----------------------------------------------------------------------------------------------
-- Calculating the Lyubeznik numbers-----
---ERROR: Still not computing values when i=j correctly except for highest Lyubeznik #
lyubeznikNumber = method()
lyubeznikNumber( ZZ, ZZ, Ideal, Ring ) := ZZ => ( i, j, I, R ) ->
(
    n := dim R;
    d := dim( R/I );
    p := char R;
    if i < 0 or j < 0 or j < i or i > d or j > d then return 0;
--    if i == 0 then i = 2;
    m := ideal R_*;
    LC := localCohomology( n-j, I, R );
    r := root LC;
    if r == 0 then return 0;
    frs1 := randomFilterRegSeq( i+1, m, r ); --frs1 has one more elm than frs
    frs := drop( frs1, -1 );
    c1 := lowerLimit(frs1,R);
    c := lowerLimit(frs,R);    
    K := matrix entries relations r;
    F := target K;
    g1 := last frs1;
    g := if #frs==0 then 0_R else last frs;
    P := ( c1*F + image K ) : g1;
    Q := c*F + (g)*F + image K;
    P1 := ((frobenius c1)*F + image frobenius K) : (g1)^p;
    Q1 := (frobenius c)*F + (g)^p*F + image frobenius K;
    M := P/Q; --ker inducedMap( F/P, F/Q );
    FM := P1/Q1; -- ker inducedMap( F/P1, F/Q1 );
    pii := product( frs, x -> x^(p-1) );
    U := matrix entries LC#generatingMorphism;
    N := makeFModule inducedMap( FM, M, pii*U );
    root N;
    degree Hom( R^1/m, root N )
)

lyubeznikTable = method()

--ADJUSTING THIS TO ONLY DO ELMS NOT ON THE DIAGONAL, just comment out the line above the comment
--and un-comment-out the comment
lyubeznikTable ( Ideal, Ring ) := Matrix => ( I, R ) ->
(
    d := dim( R/I );
    LT := apply( toList(0..d), i ->  
        apply( toList(0..d), j ->
            try lyubeznikNumber( i, j, I, R ) else -1
        )
    );
    matrix LT
)
