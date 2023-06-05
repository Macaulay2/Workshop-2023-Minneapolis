newPackage(
          "Brackets",
          Version => "0.1",
          Date => "May 5, 2023",
          Headline => "Brackets, Grassmann-Cayley Algebra, and Projective Geometry",
          Authors => {{ Name => "Tim Duff", Email => "timduff@uw.edu", HomePage => "https://timduff35.github.io/timduff35/"}},
	  PackageImports => {"SubalgebraBases"},
          AuxiliaryFiles => false,
          DebuggingMode => false
          )

export {"bracketRing", "BracketRing", "GCAlgebra", "normalForm", "gc"}

-* Code section *-

-- easy inputting of brackets and other helper functions
ZZ ZZ := (a, b) -> [a, b]
Array ZZ := (a, b) -> a | [b]
ZZ Array := (a, b) -> [a] | b
RingElement RingElement := (a, b) -> [a, b]
Array RingElement := (a, b) -> a | [b]
RingElement Array := (a, b) -> [a] | b
increment = L -> L/(l-> l+1)
sgn = sigma -> ( -- sign of a permutation
    n := length sigma;
    product(0..n-1, i -> product(i..n-1, j -> if sigma#i > sigma#j then (-1) else 1))
    )

-- AbstractGCRing is a parent class for BracketRing and GCAlgebra
AbstractGCRing = new Type of HashTable
net AbstractGCRing := G -> error "not implemented"
ring AbstractGCRing := G -> G#ring
use AbstractGCRing := G -> use ring G

-- class declaration for BracketRing
BracketRing = new Type of AbstractGCRing
-- constructor
bracketRing = method(Options => {Strategy => GroebnerBasis})
bracketRing AbstractGCRing := G -> error "not implemented"
bracketRing (List, ZZ) := o -> (vectorSymbols, d) -> (
    n := length vectorSymbols;
    if not (n >= d) then error("The first argument n in bracketRing(n, d) (representing the number of rows) is assumed to be at least the second argument d (representing the number of columns)");
    x := symbol x;
    R := QQ[x_(1,1)..x_(n,d)];
    X := matrix for i from 1 to n list for j from 1 to d list x_(i,j);
    n'choose'd := rsort(sort \ subsets(vectorSymbols, d)); -- important for "Tableux order"
    n'choose'd'Indices := rsort(sort \ subsets(#vectorSymbols, d));
    minorsX := apply(n'choose'd'Indices, R -> det X^R);
    y := symbol y; 
    bracketVariables := apply(n'choose'd, S -> y_("["|fold(S, (a, b) -> toString(a)|toString(b))|"]"));
    S := QQ[gens R, bracketVariables, MonomialOrder => {Eliminate(numgens R), GRevLex}]; -- important for "Tableaux order"
    lookupTable := new HashTable from apply(binomial(n, d), i -> increment n'choose'd'Indices#i => (gens S)#(numgens R+i));
    I := ideal apply(minorsX, bracketVariables, (m, b) -> sub(m, S) - b_S);
    ret := new BracketRing from {numrows => n, numcols => d, ring => S, ideal => I, table => lookupTable, cache => new CacheTable from {}};
    if o#Strategy === GroebnerBasis then (
	-- TODO: it's likely more efficient to apply forceGB to a known Groebner basis (Pluecker relations? van der Waerden syzygies?)
	-- TODO: allow computing with SubalgebraBases instead of Groebner bases
	G := groebnerBasis I;
	ret.cache#gb = G;
	ret.cache#syz = selectInSubring(1, G);
	) 
    else if o#Strategy === Grassmannian then (
	-- use the function "Grassmannian" to simplify construction of R, I, etc, above
	error "not implemented";
	)
    else if o#Strategy === "vanDerWaerden" then (
	error "not implemented";
	)
    else if o#Strategy === sagbi then (
	-- use SubalgebraBases package
	error "not implemented";
	)
    else error "Strategy option not recognized";
    ret
    )
bracketRing (ZZ, ZZ) := o -> (n, d) -> bracketRing(toList(1..n), d, o)
-- printing for BracketRing
net BracketRing := B -> net((symbol B)_(B#numcols, B#numrows))
-- getters for BracketRing
numrows BracketRing := B -> B#numrows 
numcols BracketRing := B -> B#numcols 
ideal BracketRing := B -> B#ideal
bracketRing BracketRing := o -> B -> B

matrix BracketRing := o-> B -> transpose genericMatrix(ring B,3,6)


-- class declaration for GCAlgebra
GCAlgebra = new Type of AbstractGCRing
gc = method(Options => {Strategy => GroebnerBasis})
-- constructor
gc (VisibleList, ZZ) := o -> (vectorSymbols, d) -> (
    n := # vectorSymbols;
    (inputMode, vectorVariables) := (
	if all(vectorSymbols, s -> instance(s, IndexedVariable) or instance(s, Symbol)) then (ZZ, vectorSymbols)
	else if all(vectorSymbols, s -> instance(s, ZZ)) and n == # unique vectorSymbols then (a := symbol a; (RingElement, apply(vectorSymbols, s -> a_s)))
	else error "incorrect input"
	);
    Bnd := bracketRing(toList vectorSymbols, d, o);
    n'Choose'd'plus1 := subsets(n,d+1);
    R := (ring Bnd)[vectorVariables, SkewCommutative=>true];
    S := R/(ideal apply(n'Choose'd'plus1, S -> product(S, s -> (vectorVariables#s)_R)));
    new GCAlgebra from {bracketRing => Bnd, ring => S}
    )

bracketRing GCAlgebra := o -> Gnd -> Gnd#bracketRing
gens GCAlgebra := Gnd -> gens ring Gnd
net GCAlgebra := Gnd -> (
    R := ring Gnd;
    n := numgens R;
    "Grassmann-Cayley Algebra generated by 1-extensors " | toString(R_0) | ".." | toString(R_(n-1))
    )
numgens GCAlgebra := Gnd -> numgens ring Gnd


GCExpression = new Type of HashTable
Bracket = new Type of GCExpression
-- pretty printing
net GCExpression := b -> (
    r := b#RingElement;
    B := bracketRing b;
    xs := drop(gens ring r, numcols B * numrows B);
    rStr := toString r;
    replace("y_","",rStr)
    )
Array _ AbstractGCRing := (A, R) -> (
    assert(#A == 1); -- For now, this function assumes the first argument A is a single-element, doubly-nested array of the form [[i1 i2 ... id]]
    A0 := A#0;
    A1 := sort toList A0;
    B := bracketRing R;
    assert(
	-- checking that the input is valid
	(instance(R, GCAlgebra) and all(A1, x -> instance(x, ambient ring R))) or 
	all(A1, i -> (instance(i, ZZ) and i >= 1 and i <= B#numrows and #A1 == B#numcols))
	);
    rowSet := first select(1, keys B#table, k -> toList A1 == k);
    new Bracket from {RingElement => (sgn A0) * B#table#rowSet, ring => R}
)


bracketRing GCExpression := o -> b -> bracketRing ring b
commonRing (GCExpression, GCExpression) := (b1, b2) -> (
    (G1, G2) := (ring b1, ring b2);
    if (G1 === G2 or G2 === bracketRing G1) then G1 else if G1 === bracketRing G2 then G2 error "Common abstract GC ring not found"
    )
degree GCExpression := A -> degree A#RingElement
ring GCExpression := b -> b#ring

RingElement _ AbstractGCRing := (b, R) -> new GCExpression from {RingElement => b, ring => R}

-- piggybacking on operators for the associated RingElement
GCExpression + GCExpression := (b1, b2) -> (
    R := commonRing(b1, b2);
    b := b1#RingElement + b2#RingElement;
    b_R
    )
GCExpression * GCExpression := (b1, b2) -> (
    R := commonRing(b1, b2);
    b := b1#RingElement * b2#RingElement;
    bR := b_R;
    if (instance(R, GCAlgebra) and (isBottomDegree bR or isTopDegree bR) and isExtensor bR) then bR_(bracketRing R) else bR
    )
GCExpression - GCExpression := (b1, b2) -> (
    R := commonRing(b1, b2);
    b := b1#RingElement - b2#RingElement;
    b_R
    )
GCExpression ^ ZZ := (b, k) -> new GCExpression from {RingElement => (b#RingElement)^k, ring => ring b}
RingElement * GCExpression := (c, b) -> new GCExpression from {RingElement => c * b#RingElement, ring => ring b}
GCExpression * RingElement := (b, c) -> new GCExpression from {RingElement => c * b#RingElement, ring => ring b}
Number * GCExpression := (c, b) -> new GCExpression from {RingElement => c * b#RingElement, ring => ring b}
GCExpression * Number := (b, c) -> new GCExpression from {RingElement => c * b#RingElement, ring => ring b}

isExtensor = method()
isExtensor GCExpression := A -> (
    assert(instance(ring A, GCAlgebra));
    1 == length terms someTerms(A#RingElement, 0, 2)
    )
isTopDegree = method()
isTopDegree GCExpression := A -> (
    assert(instance(ring A, GCAlgebra));
    (k, d) := (first degree A, numcols bracketRing A);
    (isHomogeneous A#RingElement and k == d)
    )
isBottomDegree = method()
isBottomDegree GCExpression := A -> (
    assert(instance(ring A, GCAlgebra));
    (k, d) := (first degree A, numcols bracketRing A);
    (isHomogeneous A#RingElement and k == 0)
    )
extensorSupportIndices = method()
extensorSupportIndices GCExpression := A -> (
    assert(isExtensor A);
    increment positions(first exponents A#RingElement, p -> p == 1)
    );
extensorToBracket = method()
extensorToBracket GCExpression := A -> (
    assert(instance(ring A, GCAlgebra));
    elemA := A#RingElement;
    Bnd := bracketRing A;
    local ret;
    if elemA == 0 then ret = 0 else (
	assert(first degree elemA == numcols bracketRing A);
	S := extensorSupportIndices A;
	ret = (leadCoefficient elemA) * sub(Bnd#table#S, ring elemA)
	);
    ret
    );

shuffleProduct = (A, B) -> (
    assert(isExtensor A and isExtensor B);
    Bnd := bracketRing A;
    (j, k, d) := (first degree A, first degree B, numcols bracketRing Bnd);
    if not (j+k >= d) then error "shuffle product undefined";
    P1 := permutations toList(0..d-k-1);
    P2 := permutations toList(d-k..j-1);
    shuffles := flatten flatten apply(P1, p1 -> apply(P2, p2 -> {p1 | p2, p2 | p1}));
    elemA := A#RingElement;
    suppA := support elemA;
    sum(shuffles, sigma -> (
	    -- compute the result of the shuffle product
	    bracketCoeff := extensorToBracket((leadCoefficient A#RingElement) * product(0..d-k-1, i -> suppA#(sigma#i)) * B);
	    mon := product(d-k..j-1, i -> suppA#(sigma#i));
	    ret := (sgn sigma) * bracketCoeff * mon;
	    ret
	    )
	)
    )

GCExpression _ BracketRing := (b, B) -> (
    assert(B === bracketRing b);
    bBracket := if (isTopDegree b) then lift(extensorToBracket b, ring B) else if (isBottomDegree b) then lift(b#RingElement, ring B) else error "must be an extensor of step 0 or d";
    bBracket_B
    )

GCExpression ^ GCExpression := (f, g) -> (
    G := ring f;
    assert(instance(G, GCAlgebra) and G === ring g);
    B := bracketRing G;
    fExtensors := apply(terms f#RingElement, t -> t_G);
    gExtensors := apply(terms g#RingElement, t -> t_G);
    result := (sum flatten apply(fExtensors, fE -> apply(gExtensors, gE -> shuffleProduct(fE, gE))))_G;
    if ((isBottomDegree result or isTopDegree result) and isExtensor result) then result = result_B;
    result
)
GCExpression ^ RingElement := (f, g) -> f ^ (g_(ring f))
RingElement ^ GCExpression := (f, g) -> (f_(ring g)) ^ g

normalForm = method()
normalForm GCExpression := b -> (
    assert(instance(ring b, BracketRing)); -- should we eventually allow for Grassmann Cayley algebra elements here too?
    B := bracketRing b;
    I := ideal bracketRing b;
    r := b#RingElement;
    nf := (r % I);
    nf_B
    )    

factor GCExpression := o -> g -> (
    F := factor g#RingElement;
    G := ring g;
    apply(toList F, fi -> (fi#0^(fi#1))_G)
    )


-* Documentation section *-
beginDocumentation()

doc ///
Key
  Brackets
Headline
  Brackets, Grassmann-Cayley Algebra, and Projective Geometry
Description
  Text
    Todo: add a description!
  Example
    1+1
Acknowledgement
  We thank all project contributors, the organizers of the 2023 Macaulay2 workshop in Minneapolis, IMA staff, and acknowledge support from the National Science Foundation grant DMS 2302476.
References
  Sturmfels, Bernd. <i>Algorithms in invariant theory</i>. Springer Science & Business Media, 2008.
///

-* Test section *-
TEST /// -* Sturmfels Example 3.1.10 *-
B = bracketRing(6, 3)
T = [1 4 5]_B * [1 5 6]_B * [2 3 4]_B
n = normalForm T 
assert(net n == "[256]*[145]*[134]-[356]*[145]*[124]+[456]*[145]*[123]")
-- Note: this is the same normal form Sturmfels gives, but monomials and their exponents are reverse-sorted
///

TEST /// -* Sturmfels Example 3.3.3 *-
G = gc(a..f, 3)
A = (a * d)_G
B = (b * e)_G
AB = A ^ B 
C = (c * f)_G
D = AB ^ C -- Output "2*[bde]*[acf]-2*[cdf]*[abe]" is consistent with the book's answer up to sorting and sign.
assert(net D == "2*[bde]*[acf]-2*[cdf]*[abe]")
///

TEST /// -* Desargues' Theorem*-
restart
loadPackage "Brackets"
G = gc(a..f,3)
abLine = (a * b)_G
deLine = (d * e)_G
bcLine = (b * c)_G
efLine = (e * f)_G
acLine = (a * c)_G
dfLine = (d * f)_G
pt1 = abLine ^ deLine
pt2 = bcLine ^ efLine
pt3 = acLine ^ dfLine
linePerspective = pt1 * pt2 * pt3
adLine = (a * d)_G
beLine = (b * e)_G
cfLine = (c * f)_G
pointPerspective =  adLine ^ beLine ^ cfLine
assert(net pointPerspective == "2*[bde]*[acf]-2*[cdf]*[abe]")
(n1, n2) = (normalForm pointPerspective, normalForm linePerspective);
(f1, f2) = (factor n1, factor n2)
assert(net f1#0 == "[bdf]*[ace]-[bef]*[acd]-[cdf]*[abe]-[def]*[abc]")
assert(net f2#2 == "[bdf]*[ace]-[bef]*[acd]-[cdf]*[abe]-[def]*[abc]")
///

end--

-* Development section *-
restart
loadPackage "Brackets"
check "Brackets"

uninstallPackage "Brackets"
restart
installPackage "Brackets"
viewHelp "Brackets"

