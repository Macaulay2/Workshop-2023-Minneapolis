newPackage(
          "Brackets",
          Version => "0.1",
          Date => "June 7, 2023",
          Headline => "Brackets, Grassmann-Cayley Algebra, and Projective Geometry",
          Authors => {
	      { Name => "Dalton Bidleman", Email => "", HomePage => ""},
	      { Name => "Tim Duff", Email => "timduff@uw.edu", HomePage => "https://timduff35.github.io/timduff35/"},
	      { Name => "Jack Kendrick", Email => "", HomePage => ""},
	      { Name => "Thomas Yahl", Email => "", HomePage => ""},
	      { Name => "Michael Zeng", Email => "", HomePage => ""}		      
	      },
	  PackageImports => {"SubalgebraBases"},
          AuxiliaryFiles => false,
          DebuggingMode => true
          )

export {"AbstractGCRing", "bracketRing", "BracketRing", "GCAlgebra", "normalForm", "gc", "toBracketPolynomial"}

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
AbstractGCRing Array := (G, A) -> (
    R := ring G;
    new AbstractGCRing from {ring => R A, cache => new CacheTable from {}}
    )


-- class declaration for BracketRing
BracketRing = new Type of AbstractGCRing
-- constructor
bracketRing = method(Options => {Strategy => GroebnerBasis})
bracketRing AbstractGCRing := G -> error "not implemented"
bracketRing (VisibleList, ZZ) := o -> (vectorSymbols, d) -> (
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
	-- TODO: it's likely more efficient to apply forceGB to a known Groebner basis (Plücker relations? van der Waerden syzygies?)
	-- TODO: allow computing with SubalgebraBases instead of Groebner bases
-*
	G := groebnerBasis I;
	ret.cache#gb = G;
	ret.cache#syz = selectInSubring(1, G);
*-
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
ZZ _ AbstractGCRing := (k, G) -> sub(k, ring G)
matrix BracketRing := o -> B -> transpose genericMatrix(ring B,numcols B, numrows B)

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
gens GCAlgebra := o -> Gnd -> gens ring Gnd
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
    B := bracketRing R;
    if not(B#numcols == #A0) then error "not enough symbols in bracket";
    A1 := sort toList A0;
    if any(A1, s -> instance(s, Symbol)) then A1 = apply(A1, g -> (g_R)#RingElement);
    assert(
	-- checking that the input is valid
	-- TODO: include better error messages
	(instance(R, GCAlgebra) and all(A1, x -> instance(x, ambient ring R))) or 
	all(A1, i -> (instance(i, ZZ) and i >= 1 and i <= B#numrows))
	);
    rowSet := first select(1, keys B#table, k -> toList A1 == k);
    new Bracket from {RingElement => (sgn A0) * B#table#rowSet, ring => R}
)


bracketRing GCExpression := o -> b -> bracketRing ring b
commonRing (GCExpression, GCExpression) := (b1, b2) -> (
    -- returns: either GCAlgebra or BracketRing in which b1 & b2 both make s
    (G1, G2) := (ring b1, ring b2);
    if (G1 === G2 or G2 === bracketRing G1) then G1 else if G1 === bracketRing G2 then G2 error "Common abstract GC ring not found"
    )
degree GCExpression := A -> degree A#RingElement
ring GCExpression := b -> b#ring

RingElement _ AbstractGCRing := (b, R) -> new GCExpression from {RingElement => b, ring => R}

-- piggybacking on operators for the associated RingElement
terms GCExpression := b -> (bRTerms := terms(b#RingElement);
     for term in bRTerms list term_(ring b))


GCExpression + GCExpression := (b1, b2) -> (
    R := commonRing(b1, b2);
    b := b1#RingElement + b2#RingElement;
    b_R
    )
GCExpression * GCExpression := (b1, b2) -> (
    R := commonRing(b1, b2);
    b := b1#RingElement * b2#RingElement;
    bR := b_R;
    bRTerms := terms bR;
    if (instance(R, GCAlgebra) and (all(bRTerms, isTopDegree) or all(bRTerms, isBottomDegree))) then sum(bRTerms, t -> t_(bracketRing R)) else bR
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

toBracketPolynomial = method();
toBracketPolynomial(RingElement, BracketRing) := (f, G) -> ( --input: polynomial, bracketring
    I := G#ideal;
    (f % I) _ G
)


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
	assert(first degree elemA == numcols bracketRing A); -- isTopDegree?
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

matrix (AbstractGCRing, List) := o -> (G, L) -> (
    listM := for i from 0 to (#L-1) list(apply(L#i, t-> t#RingElement));
    matrix listM)

-* Documentation section *-
beginDocumentation()

doc ///
Key
  Brackets
Headline
  Brackets, Grassmann-Cayley Algebra, and Projective Geometry
Description
  Text
      Fix integers $n \ge d \ge 1,$ and let $X = (x_{i j})$ be an $n\times d$ matrix of distinct variables in the polynomial ring $k [x_{i j}]$ over a fixed field $k.$    
      Think of each row of $X$ is a point in the projective space $\mathbb{P}^{d-1}$ of dimension $(d-1)$ over $k,$ so that $X$ as represents a configuration of $n$ points in this projective space.
      Many interesting geometric properties of this point configuration can be expressed in terms of the maximal minors of $X.$
      
      For notational convenience, it is common to write these minors in bracket notation.
      A bracket is an expression $[\lambda_1 \lambda_2 \ldots \lambda_d]$ representing the minor of $X$ whose rows are given by indices $1\le \lambda_1 < \lambda_2 < \ldots < \lambda_d \le n.$
    
      Formally, we may consider the map of polynomial rings
      $$\psi_{n,d} : k \left[ [\lambda_{i_1} \cdots \lambda_{i_d}] \mid 1 \le i_1 < \ldots < i_d \le n \right] \to k [X],$$
      $$ [\lambda_{i_1} \cdots \lambda_{i_d}] \mapsto \det \begin{pmatrix} x_{i_1, 1} & \cdots & x_{i_1, d} \\ \vdots & & \vdots & \\ x_{i_d 1} & \cdots & x_{i_d d}\end{pmatrix}. $$
      
      The classical bracket ring $B_{n,d}$ is the image of this map.
      This is the homogeneous coordinate ring of the Grassmannian of $(n-1)$-dimensional planes in $\mathbb{P}^{d-1}$ under its Plücker embedding.
Acknowledgement
  We thank all project contributors, the organizers of the 2023 Macaulay2 workshop in Minneapolis, IMA staff, and acknowledge support from the National Science Foundation grant DMS 2302476.
References
  Sturmfels, Bernd. {\it Algorithms in invariant theory}. Springer Science & Business Media, 2008.
///

doc ///
Key
  BracketRing
Description
  Text
    An object of class BracketRing represents the bracket ring $B_{n,d}$.
    For example, let $n=6, d=2,$ so that
      $$X=\begin{pmatrix}
        x_{1,1}&x_{1,2}\\
        x_{2,1}&x_{2,2}\\
        x_{3,1}&x_{3,2}\\
        x_{4,1}&x_{4,2}\\
        x_{5,1}&x_{5,2}\\
        x_{6,1}&x_{6,2}
        \end{pmatrix}.$$
      There are $6=\binom{4}{2}$ brackets, and the matrix $X$ represents a configuration of $6$ points on the projective line $\mathbb{P}^1.$
      These brackets are not algebraically independent, as they satisfy the quadratic Plücker relation,
      $$
      [1 2] [3 4] - [1 3] [2 4] + [1 4] [2 3] = 0.
      $$
      Some basic syntax for working with objects of class BracketRing is illustrated in the documentation page @TO bracketRing@.
///

doc ///
Key
  bracketRing
Headline
  Constructor for bracket rings
Usage
  B = bracketRing(n, d)
  B = bracketRing(vectorSymbols, d)
  B = bracketRing B'
  B = bracketRing G
Inputs
  B':BracketRing
  G:GCAlgebra
  vectorSymbols:List
  n:ZZ
  d:ZZ
Outputs
  B:BracketRing
Description
  Text
    To construct the bracket ring $B_{n,d}$ it is enough to specify two integers.
    In that case, each bracket will contain $n$ symbols which are integers between $1$ and $d.$
  Example
    B = bracketRing(6, 3)
    T = [1 4 5]_B * [1 5 6]_B * [2 3 4]_B
  Text
    One may also provide @ ofClass{VisibleList} @ of $n$ symbols.
  Example
    B2 = bracketRing(a..f, 3)
  Text
    Additionally, brackets can be interpreted as the top-degree elements of the Grassmann-Cayley algebra.
  Example
    G = gc(a..f, 3)
    B3 = bracketRing G
  Text
    See also @TO BracketRing@.
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
needsPackage "Brackets"
check "Brackets"

uninstallPackage "Brackets"
restart
installPackage("Brackets", RemakeAllDocumentation => true)
needsPackage "Brackets"
viewHelp "Brackets"



restart
needsPackage "Brackets"
B = bracketRing(6,3)
r = [1 2 3]_B
s = 2 * r
end

G = gc(a..f,3)
Glu = G [l, u]

restart
needsPackage "Brackets"
B = bracketRing(8,4)
