newPackage(
        "ASMVarieties",
        Version => "0.1", 
        Date => "",
        Authors => {
	    {Name => "Ayah Almousa", 
                Email => "almou007@umn.edu", 
                HomePage => "http://sites.google.com/view/ayah-almousa"},
            {Name => "Patricia Klein", 
                Email => "pjklein@tamu.edu", 
                HomePage => " "}
        },
        Headline => "functions for investigating ASM and matrix Schubert varieties",
        PackageExports => {
            "Depth",
	    "SimplicialComplexes",
	    "SimplicialDecomposability",
	    "Posets"
            },
        DebuggingMode => true
        )

export{
    "isASM",
    "essBoxes",
    "fultonGens",
    "schubertDetIdeal",
    "diagInit",
    "antiDiagInit"
    }

-- Utility routines --

--------------------------------
--auxiliary function for generating a generic matrix of z variables
--INPUT: integers n,m
--OUTPUT: an n by m generic matrix eith entries z_(i,j)
--NOTE: the ring automatically comes equipped with the antidiagonal term order
--TODO: allow user to input the field they want as an option
-----------------------------------
genMat := (n,m) -> (
    k := QQ;
    zEntries := flatten table (n,m,(i,j) -> (i,j));
    z := local z;
    degs := apply(zEntries,i-> i_1-i_0 + m); --are there better ways to make the antidiagonal weights? prob
    Q:=k[z_(1,1)..z_(n,m),MonomialOrder=>Weights=>degs];
    Mmut:=mutableMatrix(Q,n,m);
    for box in zEntries do (
        Mmut_(box) = z_(box_0+1, box_1+1);
        );
    matrix Mmut
)


findIndices = method()
findIndices(List, List) := (L, Bs) -> for ell in L list position(Bs, b -> b == ell)

---------------------------------------
-- Part 1. Constructing ASM Varieties
---------------------------------------

-------------------------------------
--checks if a given matrix is an ASM
--INPUT: matrix A
--OUTPUT: true if A is an ASM, false otherwise
-------------------------------------
isASM = method()
isASM Matrix := Boolean => (A) -> (
    n := numrows A;
    AAt :=  A |(transpose matrix{toList(n:-1)})| (transpose A) | (transpose matrix{toList(n:-1)});
    L0 := flatten entries AAt;
    L := delete(0, L0);
    partialSums := for i from 0 to length L-1 list(sum L_{0..i});
    (unique partialSums) == {1,0}
)

-----------------------
--INPUT: a list w corresponding to a permutation in 1-line notation 
--OUTPUT: a list of essential boxes
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
-----------------------
essBoxes = method()
essBoxes List := List => (w) -> (
    Zentries := set flatten table(#w, #w, (i,j)->(i+1,j+1)); --table of all boxes
    ones := apply(w, i->(i,position(w,j -> j==i)+1)); --locations of 1's in perm matrix
    L := new MutableList;
    for i from 0 to #w-1 do(
        for j from ones_i_0 to (#w) do (
            L#(#L) = (j, ones_i_1); --deathrays going down
            );
        for l from ones_i_1+1 to (#w) do (
            L#(#L) = (ones_i_0,l); --deathrays to the right
            );
        );
    toList (Zentries - set unique toList L) --determine survivors of death rays
    )

--TODO: add tests for essBoxes

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: list of fulton generators for matrix determinantal ideal w
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
----------------------------------------
fultonGens = method()
fultonGens List := List => (w) -> (
    zMatrix := genMat(#w,#w); --generic matrix
    ones := apply(w, i->(i,position(w,j -> j==i)+1)); --locations of 1's in perm matrix
    zBoxes := apply(essBoxes(w), i-> flatten table(i_0,i_1, (j,k)->(j+1,k+1))); --smaller matrix indices for each essential box
    ranks := apply(essBoxes(w), i-> #((set(zBoxes_(position(essBoxes(w), j-> j==i)))*(set ones)))); --ranks for each essential box
    fultonGens := new MutableList;
    for box in essBoxes(w) do (
        pos := position(essBoxes(w), i-> i==box);
        fultonGens#(#fultonGens) = (minors(ranks_pos+1, zMatrix^{0..(box_0-1)}_{0..(box_1-1)}))_*;
        );
    unique flatten toList fultonGens
    );

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: Schubert determinantal ideal for w
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
--WARNING: if you use the identity permutation your ring will be ZZ instead of Q and idk how to fix this
----------------------------------------
schubertDetIdeal = method()
schubertDetIdeal List := Ideal => (w) -> (
    ideal fultonGens w
    );

----------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: single-
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
--TODO: rename variables?
--WARNING: if you use the identity permutation your ring will be ZZ instead of Q and idk how to fix this
----------------------------
grothendieckPoly = method()
grothendieckPoly(List) := (w) -> (
    I := schubertDetIdeal w;
    Q := newRing(ring I, DegreeRank=> #w);
    numerator reduceHilbert hilbertSeries sub(I,Q)
    )

--TODO: add tests
--TODO: double grothendieck. Problem: can't rename variables

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: ANTIdiagonal initial ideal of Schubert determinantal ideal for w
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
--WARNING: This method does not like the identity permutation
----------------------------------------
antiDiagInit = method()
antiDiagInit List := monomialIdeal => (w) -> (
    monomialIdeal leadTerm schubertDetIdeal w
    );

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: diagonal initial ideal of Schubert determinantal ideal for w
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
--WARNING: This method does not like the identity permutation
----------------------------------------
diagInit = method()
diagInit List := monomialIdeal => (w) -> (
    I:= schubertDetIdeal w;
    R:= newRing(ring I, MonomialOrder=>Lex); --making new ring with diagonal term order (lex suffices)
    monomialIdeal leadTerm sub(I,R)
    )

-------------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: subword complex associated to w (i.e. SR-ideal of antiDiagInit)
--TODO: extend to more general subword complexes from Coxeter groups, not just these?
-------------------------------------------
subwordComplex = method()
subwordComplex List := simplicialComplex => (w) -> (
    simplicialComplex antiDiagInit w
    );

----------------------------------------
-- Part 2. Invariants of ASM Varieties
----------------------------------------


---------------------------------
---------------------------------
-- **DOCUMENTATION SECTION** --
---------------------------------
---------------------------------

beginDocumentation()

doc ///
    Key
      ASMVarieties
    Headline
      a package for investigating ASM Varieties
    Description
      Text
       This package provides functions for constructing and investigating ASM varieties. 
      Example
       isASM matrix{{1,0},{0,1}}
      Text
      	  In addition, this package provides functions for studying homological invariants of ASM varieties.
      Example
      	  2+3
///

doc ///
    Key
    	isASM
	(isASM, Matrix)
    Headline
    	Checks if a matrix is an alternating sign matrix.
    Usage
    	isASM(M)
    Inputs
    	M:Matrix
    Description
    	Text
	 Given an integer matrix, checks that the matrix is an ASM.
	Example
	 M = matrix{{1,0},{0,1}};
	 isASM M
	    
///


doc ///
    Key
    	schubertDetIdeal
	(schubertDetIdeal, List)
    Headline
    	Computes Schubert determinantal ideal for a given permutation.
    Usage
	schubertDetIdeal(L)
    Inputs
    	M:Matrix
	L:List
    Description
    	Text
	 Given a permutation or permutation matrix, outputs the Schubert determinantal ideal
	 associated to that matrix.
	Example
	 schubertDetIdeal({1,3,2})

///
-------------------------
-------------------------
--**TESTS SECTIONS**--
-------------------------
-------------------------

--detailed test for isASM function
TEST ///
-*
  restart
  needsPackage "ASMVarieties"
*-
  
--all should be true
L = {
    matrix{{1}},
    matrix{{1,0},{0,1}},
    matrix{{0,1},{1,0}},
    matrix{{0,1,0},{0,0,1},{1,0,0}},
    matrix{{0,1,0},{1,-1,1},{0,1,0}},
    matrix{{0,1,0,0},{0,0,1,0},{1,0,0,0},{0,0,0,1}},
    matrix{{1,0,0,0},{0,0,1,0},{0,1,-1,1},{0,0,1,0}},
    matrix{{0,0,1,0},{0,1,-1,1},{1,-1,1,0},{0,1,0,0}},
    matrix{{0,0,1,0},{1,0,-1,1},{0,1,0,0},{0,0,1,0}},
    matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}},
    matrix{{0,0,1,0,0,0,0,0},{1,0,-1,0,1,0,0,0},{0,0,0,1,-1,0,0,1},{0,0,1,-1,1,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,1,-1,1,0,0,0,0},{0,0,1,0,0,0,0,0}},
    matrix{{0,0,0,0,1,0,0,0},{0,0,1,0,-1,1,0,0},{0,0,0,1,0,0,0,0},{1,0,0,-1,1,-1,1,0},{0,1,-1,1,-1,1,0,0},{0,0,0,0,1,0,0,0},{0,0,1,0,0,0,0,0},{0,0,0,0,0,0,0,1}}
    };
assert(apply(L,isASM) == toList (12:true))



T = {
    matrix{{0}},
    matrix{{-1}},
    matrix{{1,0},{0,0}},
    matrix{{0,0},{0,1}},
    matrix{{-1,1},{1,0}},
    matrix{{0,1},{1,-1}},
    matrix{{0,0,0},{0,0,0},{0,0,0}},
    matrix{{0,1,0},{1,0,1},{0,1,0}},
    matrix{{-1,1,1},{1,-1,1},{1,1,-1}},
    matrix{{0,0,1},{1,0,0},{0,1,-1}},
    matrix{{0,1,0,0},{0,0,0,1},{1,0,0,0},{0,0,1,-1}},
    matrix{{0,0,1,0},{0,0,1,0},{1,0,-1,1},{0,1,0,0}},
    matrix{{0,0,0,1},{1,0,0,0},{-1,1,1,0},{1,0,0,0}},
    matrix{{1,0,1,0,-1},{0,1,0,0,0},{0,0,0,0,1},{0,0,0,1,0},{0,0,0,0,1}},
    matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,1,-1,0},{0,0,-1,1,1}},
    matrix{{0,-1,0,1,1},{1,-1,1,-1,1},{0,1,1,0,-1},{1,1,-1,1,-1},{-1,1,0,0,1}},
    matrix{{0,0,1,0,0,0,0,0},{1,0,1,0,1,0,0,0},{0,0,0,1,-1,0,0,1},{0,0,1,-1,1,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,1,-1,1,0,0,0,0},{0,0,1,0,0,0,0,0}}
    };
assert( apply(T, isASM) == toList (17:false))
///


--test for schubertDetIdeal
--TODO: Figure out how to make this test stop failing
--TODO: make more complicated tests
TEST ///
--Example 15.4 from Miller-Sturmfels
assert(schubertDetIdeal({1,2,3}) == ideal{});
assert(schubertDetIdeal({2,1,3})) == (ideal(z_(1,1)));
assert(schubertDetIdeal({2,3,1}) == ideal {z_(1,1),z_(1,2)});
assert(schubertDetIdeal({3,2,1}) == ideal{z_(1,1),z_(1,2),z_(2,1)});
assert(schubertDetIdeal({1,3,2}) == ideal(z_(1,1)*z_(2,2)-z_(1,2)*z_(2,1)));
///

end---------------------------------------------------------------


---------------------
--Ayah's sandbox
---------------------

restart
debug needsPackage "ASMVarieties"



------------------------------------
--Development Section
------------------------------------

restart
uninstallPackage "ASMVarieties"
restart
installPackage "ASMVarieties"
restart
needsPackage "ASMVarieties"
elapsedTime check "ASMVarieties"
viewHelp "ASMVarieties"
