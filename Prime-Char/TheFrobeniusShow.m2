--********************************************************************************************
-- Functions to load
--********************************************************************************************

----------------------------------------------------------------------------------------------
-- F-Modules subgroup
----------------------------------------------------------------------------------------------

needsPackage "TestIdeals"

FrobeniusFunctor = new Type of MethodFunction

FModule = new Type of HashTable

ModuleClass = new Type of Module

Morphism = new Type of Matrix

GeneratingMorphism = new Type of Morphism

STDIsomorphism = new Type of Matrix

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

stdIsomorphism = method()

stdIsomorphism ModuleClass := STDIsomorphism => 
    ( cacheValue symbol stdIsomorphism )( M ->
(
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

stdMap Morphism := Morphism => ( cacheValue symbol stdMap )( f ->
( 
    g := stdIsomorphism source f;
    h := stdIsomorphism target f;
    morphism( inverse(h) )*morphism( f*g )
))

FF = method()

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

generatingMorphism = method()

generatingMorphism Morphism := GeneratingMorphism => f ->
(
    if f == 0 then return new GeneratingMorphism from f;
    if not isWellDefined f then 
        error "generatingMorphism: map is not well defined";
    if target f != FF( source f ) then 
        error "generatingMorphism: map does not map a module M to F(M)";
    new GeneratingMorphism from f
)

generatingMorphism Matrix := GeneratingMorphism => f -> generatingMorphism morphism f
    
makeFModule = method()

makeFModule GeneratingMorphism := FModule => g -> 
    new FModule from { generatingMorphism => g, cache => new CacheTable }

makeFModule Matrix := FModule => g -> makeFModule generatingMorphism g

localCohomologyExt := ( i, I, R ) -> 
(
    M := R^1/I;
    f := inducedMap( M, FF M );
    E := Ext^i( f, R^1 );
    makeFModule generatingMorphism E
)

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
-- HSL subgroup
----------------------------------------------------------------------------------------------


----------------------------------------------------------------------------------------------
-- Big Matrices subgroup
----------------------------------------------------------------------------------------------
needsPackage "Depth"; -- for the Gorenstein case 
needsPackage "TestIdeals";

expDecomp = (R,p,e,f) -> apply(exponents(f),
      exponent->{coefficient(R_exponent,f)*R_(exponent //p^e),exponent%p^e});
--Gets the exponent vectors of each monomial X^u of the polynomial f, and associates to u the two-element list whose
        --first entry is cX^v and second entry is w, where c is the coefficient of X^u in f and u = p^e*v + w.

simpExpDecomp = (e,f) -> (
    -- returns a MutableHashTable whose keys are (sequences representing) exponents of monomial basis vectors
    R:= ring f; 
    p:= char R; 
    dec := expDecomp(R,p,e,f);
    hashT := new MutableHashTable;
    for t in dec do (ex:=toSequence t#1; coe:=t#0;
        if hashT#?ex then 
        hashT#ex = hashT#ex+coe else hashT#ex = coe);
    hashT
    )        

-- returns number which is the sequence as a base p expansion
getIndex = (seq,p,e) -> sum for i from 0 to #seq-1 list (seq#i*p^(e*i));


-- Sparse Version
AKmatrix = (R,e,f) -> (
    p:=char R;
    d:=dim R;
    n:=p^e;
    I:=frobenius^e(ideal(vars R));
    monBasis:=(entries basis(R^1/I))#0;
    expBasis:=flatten(apply(monBasis,exponents));
    L:=for i from 0 to p^(e*d)-1 list f*(monBasis#i);
    -- OLD VERSION
    -- T#i is a hashtable whose keys are exponents & values are coeefs
    --     appearing in f*(ith monomial)
    T:=flatten (for i from 0 to p^(e*d)-1 list ( expList := simpExpDecomp(e,L#i);
       for key in keys expList list (getIndex(key,p,e), i)=>expList#key));
    map(R^(n^d),R^(n^d),T)
)

fSplittingNumberNonGor = (R,e,f) -> (
    p=char R;
    n=p^e;
    S:=R/ideal(f);
    M:=coker(sub(AKmatrix(R,e,f),S));
    I:=ideal(vars S);
    phi:=inducedMap(S^1,module(I));
    return numgens source basis(coker(Hom(M,phi)))
)



--********************************************************************************************
-- Examples 
--********************************************************************************************

----------------------------------------------------------------------------------------------
-- F-Modules subgroup
----------------------------------------------------------------------------------------------

p = 3;
r = 2;
s = 6;
R = ZZ/p[x_(1,1)..x_(r,s)]
X = matrix toList apply(1..r, i -> toList apply(1..s, j -> x_(i,j) ) )

I = minors( r, X );

cohomDim I

localCohomology( 5, I, R ) == 0


all( 0..(s-r), i -> localCohomology( i, I, R ) == 0 )

associatedPrimes localCohomology( s-r+1, I, R )

----------------------------------------------------------------------------------------------

p = 2;
n = (2,2,3);
d = #n;
variables = fold( apply(d, i -> x_(i+1,0)..x_(i+1,n#i)), (i,j) -> i|j );
S = ZZ/p[variables]
B = apply(d, i -> transpose matrix toList apply(n#i, j-> { x_(i+1,j), x_(i+1,j+1) }) );
B = fold( B, (a,b) -> a | b )

J = minors( 2, B );

cohomDim J

associatedPrimes localCohomology( 6, J, S )

----------------------------------------------------------------------------------------------
-- HSL subgroup
----------------------------------------------------------------------------------------------

restart
needsPackage "FInvariants"
needsPackage "TestIdeals"

R=ZZ/2[x1,x2,x3,x4,x5]
I=ideal((x2)^2+(x1)*(x3),(x1)*(x2)*(x4)^2+(x3)^3*(x5),(x1)^2*(x4)^2+(x2)*(x3)^2*(x5))
HSLNumber(I,3)


R=ZZ/11[x,y,z]
I=ideal(x^5+y^5+z^5)
HSLNumber(I,2)


----------------------------------------------------------------------------------------------
-- Big Matrices subgroup
----------------------------------------------------------------------------------------------
R = ZZ/2[x,y,z]; 
f = z^2-x*y;
for i from 1 to 4 do print fSplittingNumberNonGor(R,i,f);
S = ZZ/3[x,y,z]; 
g = z^2-x*y;
for i from 1 to 2 do print fSplittingNumberNonGor(S,i,g);

AKmatrix(R,1,f)

A = AKmatrix(S,3,g);
print(numRows(A), numColumns(A));
