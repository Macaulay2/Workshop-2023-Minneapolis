-- -*- coding: utf-8 -*-
newPackage(
    "Subrings",
    Version => "1.0",
    Date => "June 6, 2023",
    Authors => {
	{Name => "Casey Hill", Email => "caseybhill2@gmail.com"},
        {Name => "Trevor K. Karn", Email => "karnx018@umn.edu"},
        {Name => "Miranda Moore", Email => "moor2340@umn.edu"},
        {Name => "Christopher O'Neill", Email => "cdoneill@sdsu.edu"}},
    Headline => "a package for subrings",
    Keywords => {"Documentation"},
    DebuggingMode => true
    )


export {"Subring",
        "subring",
        "subringGenerators",
	"presentationRing",
	"presentationMap",
        "presentationIdeal",
        "toQuotientRing",
	"isSubringElement"}

Subring = new Type of HashTable


-- a method to create subrings
subring = method()
subring Matrix := genMatrix -> (
    
    -- compute presentation ring
    R := ring genMatrix;
    nGens := numgens source genMatrix;
    k := coefficientRing R;
    p := symbol p;
    P := k[p_1..p_nGens];
    
    -- compute presentation map 
    f := map(R, P, genMatrix);
    
    S := new Subring from {
	generators => genMatrix,
	ambient => R,
	
	-- presentation ring: one variable for each generator
	presentationRing => P,
      
	-- presentation map: presentation ring --> ambient ring, image(f)=S
	presentationMap => f,
	 
	cache => new CacheTable
	};
    S
    )

subring List := genList -> (
    subring matrix {genList}
    )

presentationRing = method()
presentationRing Subring := S -> (
    S#presentationRing
    )

presentationMap = method()
presentationMap Subring := S -> (
    S#presentationMap
    )

presentationIdeal = method()
presentationIdeal Subring := S -> (
    --P := presentationRing S;
    f := presentationMap S;
    return ker f; --kernel is cached automatically
    )

-- a quotient ring isomorphic to the image of the subring inside of the presentation ring
toQuotientRing = method()
toQuotientRing Subring := S -> (
    P := presentationRing S;
    I := presentationIdeal S;
    return P/I;
    )

subringGenerators = method()
subringGenerators Subring := S -> S#generators
ambient Subring := S -> S#ambient

-*
generators = method()
generators Subring := S -> S#generators
*-

-- format printing of Subring type
net Subring := S -> (
    R := ambient S;
    P := presentationRing S;
    g := flatten entries S#generators;
    genstr := "";
    if #g <= 3 then (
        genstr = toString(g_{0 .. min(2, #g-1)});
    ) else (
        genstr = "{" | toString(g_0) | ", " | toString(g_1) | ", " | toString(g_2) | ", ...}";
    );

    "Subring of " | toString(R) | " generated by " | genstr | " with presentation ring " | toString(P)
    )


-*
-membership test: (x,S) -> boolean
-equality of subrings: S_1 == S_2




--Subring == Subring := (S_1, S_2) -> (
    -- if generators S_1 in S_2 and generators S_2 in S_1 return true, else false
   
  --  g_1 := S_1#generators;
  --  g_2 := S_2#generators;
    
  --  for i from 0 to #g_1-1 do (
--	if !(isElement(g_i, S_2)) 
  --	)
     
  --  )

-- given an element of S, write it in terms of the p_i's

*-

isSubringElement = method()
isSubringElement(RingElement, Subring) := (AGOODNAME, S) -> (
    R := ambient S;
    P := presentationRing S;
    T := tensor(R, P, MonomialOrder=>Eliminate(numgens R));
    gT := vars T;
    R2T := map(T, R, gT_{0 .. numgens R - 1});
    P2T := map(T, P, gT_{numgens R .. numgens T - 1});
    ambGens := gens R;
    subringGens := subringGenerators S;
    presGens := gens P;
    graphIdealGens := for i from 0 to (-1 + numgens P) list
     (P2T(presGens_i) - R2T(subringGens_i)_0);
    I := ideal graphIdealGens;
    M := matrix {{R2T(AGOODNAME) % I}};
    selectInSubring(1, M) == M
    )


beginDocumentation()
doc ///
 Node
  Key
   Subrings
  Headline
   a package to deal with subrings
  Description
   Text
     {\em Subrings} is a package to give basic subroutines for subrings.
    Caveat
     There are other subring flavor things out there. 
    Subnodes
     Subring
     subring
     presentationRing
     presentationMap
     presentationIdeal
     toQuotientRing
    Example
     needsPackage("Subrings")
///
doc ///
    Key
        subring
        (subring, Matrix)
        (subring, List)
    Headline
         Construct a subring of a polynomial ring
    Usage
        S = subring M
        S = subring L
    Inputs
        M:Matrix
            of generators for a subring of a @ ofClass{PolynomialRing} @
        L:List
            of generators for a subring of a @ ofClass{PolynomialRing} @
    Outputs
        S:Subring
            the subring of the polynomial ring
    Description
        Text
            An easy way to specify a subring is to specify the ambient ring and the generators of a desired subring as a matrix. The ambient ring is implicit in the constructor function.
        Example
            R = QQ[x,y]
            M = matrix(R, {{x^2, x*y, y^2}})
            S = subring M
        Text
            This function also accepts a list of elements of a polynomial ring as an input.
        Example
            R = QQ[x,y]
            L = {x^2, x*y, y^2}
            S = subring L
///
doc ///
    Key
        presentationIdeal
        (presentationIdeal, Subring)
    Headline
        Compute the presentation ideal of a subring
    Usage
        I = presentationIdeal S
    Inputs
        S:Subring
            a @ ofClass{Subring} @ of a @ ofClass{PolynomialRing} @
    Outputs
        I:Ideal
            the presentation ideal of the subring
    Description
        Text
            This function finds the presentation ideal of the subring, which is defined to be the presentation ring modulo the presentation map.
        Example
            R = QQ[x,y]
            S = subring {x^2, x*y, y^2}
            I = presentationIdeal S
///

doc ///
 Node
  Key
   presentationRing
  Headline
   a polynomial ring with a variable for each subring generator
  Description
   Text
    The {\tt presentationRing} of a subring {\tt S} is a polynomial ring
    with a variable for each generator of {\tt S}.
   Example
    R = QQ[x,y];
    g = {x^2, x*y, y^2};
    S = subring g
    presentationRing S
///

doc ///
 Node
  Key
   presentationMap
  Headline
   the map from the {\tt presentationRing} to a subring {\tt S}
  Description
   Text
    There is a map sending each generator of the presentation ring to
    the corresponding element of the ambient ring. This is that map.
   Example
    R = QQ[x,y];
    g = {x^2, x*y, y^2};
    S = subring g;
    f = presentationMap S
    P = presentationRing S
    p = P_0 * P_1 - P_2
    f(p)
///


TEST ///
    assert(true == true)
///

end--

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
