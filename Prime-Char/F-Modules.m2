loadPackage "TestIdeals"

----------------------------------------------------------------------------------------------
-- Auxiliary Functions
----------------------------------------------------------------------------------------------

identityMatrix := n -> id_( ZZ^n )

----------------------------------------------------------------------------------------------
-- Types
----------------------------------------------------------------------------------------------

FrobeniusFunctor = new Type of MethodFunction

FModule = new Type of HashTable

ModuleClass = new Type of Module

Morphism = new Type of Matrix

GeneratingMorphism = new Type of Morphism

STDIsomorphism = new Type of Matrix

----------------------------------------------------------------------------------------------
-- ModuleClass and Morphism
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

morphism = method()

morphism Matrix := Morphism => f -> 
( 
    M := moduleClass source f;
    N := moduleClass target f;
    new Morphism from map( N, M, matrix f )
)

-- Putting modules and homomorphisms in a form that can be easily frobenified    

stdIsomorphism = method()

-- stdIsomorphism(M) returns an isomorphism M' --> M, where M' is a quotient of a free module 
stdIsomorphism ModuleClass := STDIsomorphism => 
    ( cacheValue symbol stdIsomorphism )( M ->
(
    -- If M is already a quotient of a free module, don't mess with it.
    -- BUT... maybe it would be good to simplify it, if possible.  
--    if isQuotientModule M then return new STDIsomorphism from id_M;    
    P := minimalPresentation M;
    local g;
    if P#cache#?pruningMap then g = P#cache#pruningMap
    else 
    (
        n := rank source generators M;
        g = map( M, P, entries identityMatrix n )
    );
    new STDIsomorphism from g 
))

stdMap = method()

-- Given a module homomophism f: M --> N, stdMap(f) returns the homomorphism f': M' --> N' 
-- induced by the standard isomorphisms M' --> M and N' --> N. That is, f' makes the following 
-- diagram commute:
--         f 
--     M  --->  N
--     ^        ^
--     |   f'   |
--     M' --->  N'

stdMap Morphism := Morphism => ( cacheValue symbol stdMap )( f ->
( 
    g := stdIsomorphism source f;
    h := stdIsomorphism target f;
    morphism( inverse(h) )*morphism( f*g )
))

----------------------------------------------------------------------------------------------
-- Frobenius Functor
----------------------------------------------------------------------------------------------

FF = method()

-- TODO: need to check if ring is a polynomial ring (or regular) 

FF ( ZZ, ModuleClass ) := ModuleClass => ( e, C ) -> 
(
    if isFreeModule C or C == 0 then return C; 
    R := ring C;
    p := char R;
    Rel := relations source stdIsomorphism C;
    -- adjustments of degrees are needed to keep things homogeneous and maps composable
    degsource := - p^e * ( degrees source Rel );
    degtarget := - p^e * ( degrees target Rel );
    new ModuleClass from coker map( R^degtarget, R^degsource, entries frobenius^e Rel )
)

FF ( ZZ, Module ) := ModuleClass => ( e, M ) -> FF( e, moduleClass M )

FF ( ZZ, Morphism ) := Morphism => ( e, g ) -> 
(
    f := stdMap g;
    new Morphism from map( FF(e, target f), FF(e, source f), frobenius^e f )
)

FF ( ZZ, Matrix ) := Morphism => ( e, f ) -> FF( e, morphism f )   

FF Thing := M -> FF( 1, M )

FF = new FrobeniusFunctor from FF

FrobeniusFunctor ^ ZZ := ( f, n ) -> ( x -> f( n, x ) )

----------------------------------------------------------------------------------------------
-- F-Modules
----------------------------------------------------------------------------------------------

generatingMorphism = method()

-- Given a standard map f, generatingMorphism(f) verifies whether f well defined, and whether
-- f maps a module M to F(M). If so, it returns a GeneratingMorphism identical to f.
generatingMorphism Morphism := GeneratingMorphism => f ->
(
    if f == 0 then return new GeneratingMorphism from f;
    if not isWellDefined f then 
        error "generatingMorphism: map is not well defined";
    if target f != FF( source f ) then 
        error "generatingMorphism: map does not map a module M to F(M)";
    new GeneratingMorphism from f
)

-- Standardizes a map and tries to make a GeneratingMorphism from it.
generatingMorphism Matrix := GeneratingMorphism => f -> generatingMorphism morphism f
    
makeFModule = method()

makeFModule GeneratingMorphism := FModule => g -> 
    new FModule from { generatingMorphism => g, cache => new CacheTable }

--- Compute a generating morphism for H_I^i(R)
localCohomology = method()

localCohomology ( ZZ, Ideal, Ring ) := FModule => ( i, I, R ) -> 
(
    -- TODO: check that I is ideal of R, positive characteristic, etc.
    M := R^1/I;
    f := inducedMap( M, FF M );
    E := Ext^i( f, R^1 );
    makeFModule generatingMorphism E
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

cohomDim = method()

cohomDim Ideal := ZZ => I ->
(
    R := ring I;
    n := #(trim I)_*;
    while localCohomology( n, I, R) == 0 do n = n-1;
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
    isFRE = isFilterRegElement(L#0,I,M,module ideal(0_R));
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
    while not isFilterRegElement(G#i,I,M,module ideal(0_R)) do (i=i+1;);
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
    if o.Homogeneous then random( deg, I, Density => density )
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

randomFilterRegSeq (ZZ,Ideal,Ring):= List => o -> (n,I,R) ->
(
    randomFilterRegSeq(n,I,R^1);
)
    


----------------------------------------------------------------------------------------------

-- generates limit closure and lower limit ideal


ascendingIdealEquality = method( Options => {MaxTries => infinity})
-----this function should eventually check to make sure the chain of ideals is ascending

-----this function specifically finds the FIRST spot that we retain equlaity and returns
-----the union of ideals up to that spot (equivalently the ideal at that spot)
ascnedingIdealEquality (Function) := Ideal => o -> f ->
(
    i:=1;
    isEqual = false;
    while (not isEqual and i<o.MaxTries) do (
	isEqual = (f(i)==f(i+1));
	i=i+1;
	);
    if not isEqual then error "ascnedingIdealEquality: Reached maximum limit of tries.";
    f(i)
)

limitClosure = method()
limitClosure(BasicList) := Ideal => L ->
(
    ascnedingIdealEquality(j -> (ideal(apply(L,x->x^(p^j)):ideal(product(apply(L,x->x^(p^j-1))))))
)

lowerLimit = method()
lowerLimit(BasicList) := Ideal => L ->
(
    LC:=limitClosure(drop(L,-1));
    ascnedingIdealEquality( j -> ideal(LC:(L#(-1)^j)))
)
    
