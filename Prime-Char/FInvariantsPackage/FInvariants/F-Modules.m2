needsPackage "TestIdeals"

needsPackage "Saturation"

----------------------------------------------------------------------------------------------
-- Auxiliary Functions
----------------------------------------------------------------------------------------------

FrobeniusFunctor = new Type of MethodFunction

FModule = new Type of HashTable

GeneratingMorphism = new Type of Matrix

----------------------------------------------------------------------------------------------
-- Frobenius Functor
----------------------------------------------------------------------------------------------

isGoodRing = method()

isGoodRing Ring := Boolean => R -> isPolynomialRing R and char R > 0
-- TODO In the future, this should check regularity instead. 

FFmethod = method()

FFmethod ( ZZ, Module ) := Module => ( e, M ) -> 
(
    if isFreeModule M or M == 0 then return M; 
    R := ring M;
    if not isGoodRing R then error "FF is only implemented for modules and morphisms over regular rings of positive characteristic";
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

FFmethod ( ZZ, Matrix ) := Matrix => ( e, f ) -> 
    map( FF(e, target f), FF(e, source f), entries frobenius^e f )

FFmethod Thing := M -> FF( 1, M )

FF = new FrobeniusFunctor from FFmethod

FrobeniusFunctor ^ ZZ := ( f, n ) -> ( x -> f( n, x ) )

----------------------------------------------------------------------------------------------
-- F-Modules
----------------------------------------------------------------------------------------------

generatingMorphism = method()

-- Given f, generatingMorphism(f) verifies whether f well defined, and whether
-- f maps a module M to F(M). If so, it returns a GeneratingMorphism identical to f.
generatingMorphism Matrix := GeneratingMorphism => f ->
(
    if f == 0 then return new GeneratingMorphism from f;
    if not isWellDefined f then 
        error "generatingMorphism: map is not well defined";
    if minimalPresentation( target f ) != minimalPresentation( FF( source f ) ) then 
        error "generatingMorphism: does not map a module M to F(M)";
    new GeneratingMorphism from f
)
    
makeFModule = method()

makeFModule GeneratingMorphism := FModule => g -> 
    new FModule from { cache => new CacheTable from { generatingMorphism => g } }

makeFModule Matrix := FModule => g -> makeFModule generatingMorphism  g

--- Compute a generating morphism for H_I^i(R)

localCohomologyExt = ( i, I ) -> 
(
    R := ring I;
    M := R^1/I;
    f := inducedMap( M, FF M );
    E := Ext^i( f, R^1 );
    makeFModule E
)

--- Compute a generating morphism for H_I^i(R)
localCohomologyFilter = ( i, I ) -> 
(
    R := ring I;
    filterSeq := randomFilterRegSeq( i, I, R );
    J := ideal filterSeq;
    p := char R;
    u := ( product filterSeq )^( p-1 );
    time K := firstEquality( 
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

localCohomology ( ZZ, Ideal ) := FModule => o -> ( i, I ) -> 
(   
    R := ring I;
    if not isGoodRing R then error "localCohomology is only implemented for regular rings of positive characteristic";
    print "oi";
    if (I#cache)#?(localCohomology, i) then return I#cache#(localCohomology, i);
    lc := if o.Strategy === Ext then localCohomologyExt( i, I )
    else localCohomologyFilter( i, I );
    I#cache#(localCohomology, i) = lc;
    lc
)

localCohomology ( ZZ, Ideal, Ring ) := FModule => o -> ( i, I, R ) ->
(
    if R =!= ring I then error "localCohomology: expected an ideal of the given ring";
    localCohomology( i, I, o )
)
  
root = method()

root FModule := Module => ( cacheValue symbol root )( M ->
(
    g := M#cache#generatingMorphism;
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
    N := (source g)/K;
    M#cache#(symbol rootMorphism) = map( FF N, N, g );  
    N
))
           
ZZ == FModule := ( n, M ) -> M == n
           
FModule == ZZ := ( M, n ) -> 
(
    if n =!= 0 then error "Attempted to compare an FModule to nonzero integer";
    -- check while generating root == 0
    ( root M )  == 0
)

cohomDim = method( Options => { Strategy => Ext } )

cohomDim Ideal := ZZ => o -> ( cacheValue symbol cohomDim )( I ->
(
    n := #(trim I)_*;
    while localCohomology( n, I, o) == 0 do n = n-1;
    n
))

associatedPrimes FModule := o -> M -> 
(
    if M#cache#?associatedPrimes then M#cache#associatedPrimes
    else
    (
        ap := associatedPrimes( root M, o );
        M#cache#associatedPrimes = ap;
        ap
    )
)

----------------------------------------------------------------------------------------------
-- Generating random generating morphisms and FModules
----------------------------------------------------------------------------------------------

randomGeneratingMorphism = method( Options=> { Degree => 1 } )

-- produces a random generating morphism M --> F(M)
randomGeneratingMorphism Module := GeneratingMorphism => o -> M ->
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

isFilterRegElement ( RingElement, Ideal, Module, Module ) := Boolean => ( x, I, M, N ) ->
(
    T := ker( x*id_( M/N ) );
    J := radical ann T;
    isSubset( I, J )
)

isFilterRegSeq = method()

isFilterRegSeq ( BasicList, Ideal, Module ) := Boolean => ( L, I, M ) ->
(
    if not isSubset( ideal L, I ) then error "isFilterRegSeq: The sequence is not contained in the ideal";
    l := #L;
    R := ring I;
    isFRE := isFilterRegElement( L#0, I, M, ( ideal 0_R )*M );
    N := ( ideal L#0 )*M;
    i := 1;
    while isFRE and i < l do 
    (
	isFRE = isFilterRegElement( L#i, I, M, N );
	N = N + ( ideal L#i )*M;
	i = i + 1;
    );
    isFRE
)

isFilterRegSeq ( BasicList, Ideal, Ring ) := Boolean => ( L, I, R ) ->
    isFilterRegSeq( L, I, module R )

isFilterRegSeq ( BasicList, Ideal, Ideal ) := Boolean => ( L, I, J ) ->
    isFilterRegSeq( L, I, module J )

filterRegSeq = method()

filterRegSeq (ZZ,Ideal,Module) := List => (n,I,M) ->
(
    --ap:= associatedPrimes M;
    --ap:= select(ap,p -> (not isSubset(I,p)))
    
    --Strategy of ordering these generators obtained from Eisenbud-Huneke-Vasconcelos.m2 in the
    ---PrimaryDecomposition package.
    R := ring I;
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
    R := ring I;
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
    randomFilterRegSeq( n, I, module R )
    
----------------------------------------------------------------------------------------------

-- generates limit closure and lower limit ideal

firstEquality = method( Options => { Tries => infinity } )

-- Given a function f defined on the positive integers, firstEquality finds the first 
-- value f(i) that is equal to its successor, f(i+1).
firstEquality Function := Ideal => o -> f ->
(
    i := 0;
    isEqual := false;
    while not isEqual and i < o.Tries do 
    (
	i = i+1;
	isEqual = f(i) == f(i+1)
    );
    if not isEqual then error "firstEquality: Reached maximum limit of tries.";
    f(i)
)

limitClosure = method()
limitClosure BasicList := Ideal => L -> 
(
    if #L == 0 then error "limitClosure: limit closure should be the 0 ideal, but cannot check the ring it lives in since the list is empty.";
    p := char ring L#0;
    firstEquality(j -> (frobenius^j(ideal(L)):ideal(product(L,x->x^(p^j-1)))))
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
    firstEquality( j -> LC : L#(-1)^j )
)

lowerLimit RingElement := Ideal => f -> lowerLimit { f }
    
----------------------------------------------------------------------------------------------
-- Calculating the Lyubeznik numbers-----
---ERROR: Still not computing values when i=j correctly except for highest Lyubeznik #
lyubeznikNumber = method()

lyubeznikNumber( ZZ, ZZ, Ideal ) := ZZ => ( i, j, I ) -> 
(
    ln := ( cacheValue ( symbol lyubeznikNumber, i, j ) )( I -> (
    R := ring I;
    n := dim R;
    d := dim( R/I );
    p := char R; 
    if i < 0 or j < 0 then error "lyubeznikNumber: expected nonnegative integers"; 
    if i > j or j > d then return 0;
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
    M := (( c1*F + image K ) : g1)/(c*F + (g)*F + image K);
    pii := product( frs, x -> x^(p-1) );
    U := matrix entries LC#cache#generatingMorphism;
    N := makeFModule inducedMap( FF M, M, pii*U );
    root N;
    degree Hom( R^1/m, root N )
    ));
    ln I
)

lyubeznikNumber ( ZZ, ZZ, Ideal, Ring ) := ( i, j, I, R ) ->
(
    if R =!= ring I then error "lyubeznikNumber: expected an ideal of the given ring";
    lyubeznikNumber( i, j, I )
)

lyubeznikTable = method()

lyubeznikTable ( Ideal, Ring ) := Matrix => ( I, R ) ->
(
    d := dim( R/I );
    LT := toList apply( 0..d, i ->  
        toList( (i:0) | apply( i..d, j -> lyubeznikNumber( i, j, I, R ) ) )
    );
    matrix LT
)


