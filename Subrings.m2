-- -*- coding: utf-8 -*-
newPackage(
    "Subrings",
    Version => "1.0",
    Date => "June 6, 2023",
    Authors => {
	{Name => "Casey Hill", Email => "casey.hill@uky.edu"},
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

-- a method to create subrings from a `Matrix` of generators
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

-- a method to create subrings from a `List` of generators
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

subringGenerators = method() -- TODO: is there a way to call just "gens"/"generators"?
subringGenerators Subring := S -> S#generators

generators Subring := Matrix => opts -> S -> (
    if  S#?generators then S#generators
    else if S.cache.?generators then S.cache.generators
    else S.cache.generators = S#generators)

ambient Subring := S -> S#ambient

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

-- given a subring S and an element x of the ambient ring, checks whether x is in S
isSubringElement = method()
isSubringElement(RingElement, Subring) := (r, S) -> (
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
    M := matrix {{R2T(r) % I}};
    selectInSubring(1, M) == M
    )

-- equality of Subrings
-- check that every generator of S1 is in S2 and every generator of S2 is in S1
Subring == Subring := (S1, S2) -> (
    if not (ambient S1) === (ambient S2) then return false;
    all(entries(subringGenerators S1)_0, f -> isSubringElement(f, S2)) and all(entries(subringGenerators S2)_0, f -> isSubringElement(f, S1))
    )



-*

-- given an element of S, write it in terms of the p_i's

*-

-----------------------------------------
--- Load documentation files ------------
-----------------------------------------
beginDocumentation()

load "./SubringsDoc.m2"
load "./SubringsTests.m2"

end--

restart
installPackage "Subrings"
R = QQ[x]
S1 = subring {x, x^2}
S2 = subring {x^2}
S3 = subring {x}
S1 == S3 --true
S1 == S2 --false

net S3

mingens S1
presentationRing S1
mingens

gens S1

-- run tests
check Subrings

