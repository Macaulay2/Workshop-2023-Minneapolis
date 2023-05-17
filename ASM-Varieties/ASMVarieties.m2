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
    "fultonGens",
    "schubertDetIdeal",
    "diagInit",
    "antiDiagInit",
    "rankMatrix",
    "essentialBoxes",
    "subwordComplex",
    "grothendieckPoly",
    "rotheDiagram"
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

--------------------------------
--auxiliary function for making a permutation matrix out of a perm in 1-line notation
--INPUT: a list w, which is a permutation in 1-line notation
--OUTPUT: a permutation matrix A corresponding to to w
--TODO: add check that w is indeed a permutation
-----------------------------------
permToMatrix := (w) -> (
    n := #w;
    (id_(ZZ^n))_(apply(w, i-> i-1))
    )

-*
findIndices(List, List) := (L, Bs) -> (for ell in L list position(Bs, b -> (b == ell)))
*-
--------------------------------------------
--------------------------------------------
--**Part 1. Constructing ASM Varieties**--
--------------------------------------------
--------------------------------------------

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

----------------------------------------
--Computes rank matrix of an ASM
--INPUT: an (n x n)- alternating sign matrix A OR a 1-line perm w
--OUTPUT: an (n x n) integer matrix of ranks of each entry
--Author: Yuyuan Luo
--TODO: add error check that A is an ASM using isASM
--TODO: add tests and documentation for this function
----------------------------------------

rankMatrix = method()
rankMatrix Matrix := Matrix => (A) -> (
    n := numrows A;
    rankA := {};
    for i from 0 to n-1 do (
	temp := toList(n:0);
        for j from 0 to n-1 do (
            if (j>0) then prev := temp#(j-1);
            if (i==0 and j==0) then temp=replace(j,A_(0,0),temp)
            else if (i==0) then temp=replace(j,prev+A_(i,j),temp)
            else if (j==0) then temp=replace(j,rankA_{i-1}#0#(j)+A_(i,j),temp)
            else temp=replace(j,rankA_{i-1}#0#(j)+prev-rankA_{i-1}#0#(j-1)+A_(i,j),temp));
        rankA=append(rankA, temp));
    matrix rankA
    )
rankMatrix List := Matrix => (w) -> (
    A := permToMatrix w;
    rankMatrix A
    )

 
-*
--Ayah's version deprecated in favor of Yuyuan's more general version
--kept for reference, we can delete later
essentialBoxes List := List => (w) -> (
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
*-

-----------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
    	--OR an alternating sign matrix A
--OUTPUT: a list of boxes in the Rothe diagram for A
--TODO: add check that list w is actually a permutation/that A is an ASM
--TODO: add documentation
--TODO: needs fixing
--Sometimes this gives the essential set and sometimes this gives all boxes
-----------------------
rotheDiagram = method()
--this version for general ASMs due to Yuyuan Luo
rotheDiagram Matrix := List => (ASM) -> (
    n := numrows ASM;
    essL := {};
    seen := {};
    for i from 0 to n-1 do (
        x := toList(n:0);
        seen = append(seen,x);
	);

    for i from 0 to n-1 do (
        for j from 0 to n-1 do (
            if (ASM_(i,j) == -1) then (essL = append(essL,{i,j}); continue;);
            if (ASM_(i,j) == 1) then (
                seenc := seen_{i}#0;
                seenc = replace(j, 1, seenc);
                seen = replace (i, seenc, seen);
                
                if (i < n-1) then (
                    for k from i+1 to n-1 do (
                        if (ASM_(k,j) == 1 or ASM_(k,j) == -1) then break;
                        seenc := seen_{k}#0;
                        seenc = replace(j,1,seenc);
                        seen = replace(k,seenc,seen);
                    );
                );
                if (j < n-1) then (
                    for k from j+1 to n-1 do (
                        if (ASM_(i,k) == 1 or ASM_(i,k) == -1) then break;
                        seenc := seen_{i}#0;
                        seenc = replace(k,1,seenc);
                        seen = replace(i,seenc,seen);
                    );
                );
                    
            );
        
        );
    );

    for i from 0 to n-1 do (
        for j from 0 to n-1 do (
            if ((seen_{i}#0#(j)) == 0) then essL = append(essL,{i,j});
        );
    );
sort unique apply(essL, i-> {i_0+1,i_1+1})
)

rotheDiagram List := List => (w) -> (
    A := permToMatrix w;
    essentialBoxes(A)
    )


-----------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
    	--OR an alternating sign matrix A
--OUTPUT: a list of essential boxes in the Rothe diagram for A
--TODO: add check that list w is actually a permutation/that A is an ASM
--TODO: add documentation
-----------------------
essentialBoxes = method()
essentialBoxes Matrix := List => (A) -> (
    boxes := rotheDiagram(A);
    badBoxes := apply(boxes, i->(positions(boxes,j->(j=={i_0,i_1+1}))|positions(boxes,j->(j=={i_0,i_1+1}))));
    essBoxes := positions(badBoxes, i-> i=={});
    boxes_essBoxes
    )
essentialBoxes List := List => (w) -> (
    essentialBoxes permToMatrix w
    )


----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: list of fulton generators for matrix determinantal ideal w
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
----------------------------------------
fultonGens = method()
fultonGens Matrix := List => (A) -> (
    zMatrix := genMat(numrows A, numcols A); --generic matrix
    rankMat := rankMatrix A; --rank matrix for A
    essBoxes := essentialBoxes A;
    zBoxes := apply(essBoxes, i-> flatten table(i_0,i_1, (j,k)->(j+1,k+1))); --smaller matrix indices for each essential box
    ranks := apply(essBoxes, i-> rankMat_(i_0-1,i_1-1)); --ranks for each essential box
    fultonGens := new MutableList;
    for box in essBoxes do (
        pos := position(essBoxes, i-> i==box);
        fultonGens#(#fultonGens) = (minors(ranks_pos+1, zMatrix^{0..(box_0-1)}_{0..(box_1-1)}))_*;
        );
    unique flatten toList fultonGens
    );
fultonGens List := List => (w) -> (
    A := permToMatrix w;
    fultonGens A
    )
----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: Schubert determinantal ideal for w
--TODO: extend to work for inputting a permutation matrix, a string w, or an ASM
--TODO: add check that list w is actually a permutation
--WARNING: if you use the identity permutation your ring will be ZZ instead of Q and idk how to fix this
----------------------------------------
schubertDetIdeal = method()
schubertDetIdeal Matrix := Ideal => (A) -> (
    ideal fultonGens A
    )
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
grothendieckPoly(Matrix):= (A) -> (
    I := schubertDetIdeal A;
    Q := newRing(ring I, DegreeRank=> numcols A);
    numerator reduceHilbert hilbertSeries sub(I,Q)
    )
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
antiDiagInit Matrix := monomialIdeal => (A) -> (
    monomialIdeal leadTerm schubertDetIdeal A
    );
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
diagInit Matrix := monomialIdeal => (A) -> (
    I:= schubertDetIdeal A;
    R:= newRing(ring I, MonomialOrder=>Lex); --making new ring with diagonal term order (lex suffices)
    monomialIdeal leadTerm sub(I,R)
    )
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
subwordComplex Matrix := simplicialComplex => (A) -> (
    simplicialComplex antiDiagInit A
    );

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
       Many of the functions in this package can take as input either a permutation matrix in 1-line notation,
       or an alternating sign matrix.
      Example
	   w = {1,5,3,4,2};
	   essentialBoxes w   
	   netList fultonGens w
      Text
       	   This package also contains functions for studying homological properties of ASM varieties.
      Example
      	  grothendieckPoly w
      	  betti res diagInit w
	  betti res antiDiagInit w	   
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
    	    M = matrix{{0,0,1,0,0,0,0,0},{1,0,1,0,1,0,0,0},{0,0,0,1,-1,0,0,1},{0,0,1,-1,1,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,1,-1,1,0,0,0,0},{0,0,1,0,0,0,0,0}}
	    isASM M
	    N = matrix{{0,-1,0,1,1},{1,-1,1,-1,1},{0,1,1,0,-1},{1,1,-1,1,-1},{-1,1,0,0,1}}
	    isASM N

	    
///


doc ///
    Key
    	schubertDetIdeal
	(schubertDetIdeal, List)
    Headline
    	Computes Schubert determinantal ideal for a given permutation.
    Usage
    	schubertDetIdeal(M)
	schubertDetIdeal(L)
    Inputs
    	M:Matrix
	L:List
    Description
    	Text
	 Given an alternating sign matrix or a permutation in 1-line notation, 
	 outputs the Schubert determinantal ideal associated to that matrix.
	Example
	 schubertDetIdeal({1,3,2})
	 schubertDetIdeal(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

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

TEST ///
--Example 2.1 in Weigandt "Prism Tableaux for ASMs"
A = matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}};
assert(isASM A)
assert(sort essentialBoxes(A) == {{1,3},{2,1},{3,2}})

///
-*
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
*-

end---------------------------------------------------------------


---------------------
--Ayah's sandbox
---------------------

restart
debug needsPackage "ASMVarieties"

M = matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}},
M = {1,4,3,5,2}
isASM M
rotheDiagram M
essentialBoxes M
grothendieckPoly M
netList fultonGens M
subwordComplex M
betti res diagInit M
betti res antiDiagInit M
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
