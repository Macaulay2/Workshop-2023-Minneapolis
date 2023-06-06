restart
debug needsPackage "Brackets"
-*

The goal of this exercise is to prove that 4 lines in P^3 can have either 1, 2, or infinitely many transversals. 

To do this, we fix 4 pairs of points in P^3, and look at the condition that the common transversal of the first three lines
also intersects with the other line.  

The first part of this script is to be converted as the Grassmann strategy for the GCExpression method. 

*-
----------------------------------------------------------------
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

----------------------------------------------------------------

I = Grassmannian(1,3)
needsPackage "Brackets"
B = bracketRing(4,2)
apply(gens ring I, v -> (
    S := last baseName v;
    A := [new Array from apply(S, i -> i + 1)];
    A_B
    )
    )




-- grassmann method for computing GCalgebra 
G37 = Grassmannian(3,7) --ideal of Pluecker relations

G37gb = gb G37  -- Groebner basis

G37ring = ring G37 -- grassmannian ring before quotient

G = G37ring/G37 -- quotient ring is desired GCalgebra
--G = gc(a .. h, 4) -- old method that is too computationally expensive

----------------------------------------------------------------
R = (G, [u,v]) -- algebra with indeterminants u and v over GCalgebra G

l1 = a*b; --
l2 = c*d;
l3 = e*f;
--l4 = g*h;

l = (((u*a + v*b)*c*d)^(e*f)) * (u*a + v*b)

quantity1 = l*g*h

toBracketPolynomial(quantity1, G)


