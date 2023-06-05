newPackage(
        "MatrixSchubert",
        Version => "0.1", 
        Date => "",
        Authors => {
	    {Name => "Ayah Almousa", 
                Email => "almou007@umn.edu", 
                HomePage => "http://sites.google.com/view/ayah-almousa"},
            {Name => "Patricia Klein", 
                Email => "pjklein@tamu.edu", 
                HomePage => " "},
	    {Name => "Yuyuan Luo",
		Email => "",
		HomePage=> ""}
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
    "isPartialASM",
    "fultonGens",
    "schubertDetIdeal",
    "diagLexInit",
    "antiDiagInit",
    "rankMatrix",
    "essentialBoxes",
    "subwordComplex",
    "grothendieckPoly",
    "rotheDiagram",
    "permToMatrix",
    "composePerms",
    "isPerm",
    "minimalRankTable"
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

--INPUT: A list w
--OUTPUT: TRUE if w is a permutation; else FALSE

isPerm = method()
isPerm List := Boolean => (w) -> (
    n := #w;
    (sort w) == toList(1..n)
    )

--------------------------------
--auxiliary function for making a permutation matrix out of a perm in 1-line notation
--INPUT: a list w, which is a permutation in 1-line notation
--OUTPUT: a permutation matrix A corresponding to to w
--NOTE: this might be the transpose of what some would want
--TODO: add documentation
-----------------------------------
permToMatrix = method()
permToMatrix List := Matrix => (w) -> (
    if not(isPerm w) then error("The input must be a permutation.");
    n := #w;
    transpose (id_(ZZ^n))_(apply(w, i-> i-1))
    )

-*
findIndices(List, List) := (L, Bs) -> (for ell in L list position(Bs, b -> (b == ell)))
*-


-----------------------------------------------
-----------------------------------------------
--**Useful functions for permutations**
-----------------------------------------------
-----------------------------------------------


------------------------------------
--INPUT: Two permutations in one-line notation
--OUTPUT: the composition of the two permutations w = vu
------------------------------------
composePerms = method()
composePerms (List, List) := List => (u,v) -> (
    if not (isPerm u) then error("The first argument is not a permutation.");
    if not (isPerm v) then error("the second argument is not a permutation.");
    u0 := apply(u, i->i-1);
    v0 := apply(v, i->i-1);
    apply(u0_v0, i-> i+1)
    )

--------------------------------------------
--------------------------------------------
--**Part 1. Constructing ASM Varieties**--
--------------------------------------------
--------------------------------------------

-------------------------------------
--checks if a given matrix is a partial ASM
--INPUT: matrix A
--OUTPUT: true if A is a partial ASM, false otherwise
-------------------------------------
isPartialASM = method()
isPartialASM Matrix := Boolean => (A) -> (
    n := numrows A;
    m := numcols A;
    rowCheck := new MutableList;
    colCheck := new MutableList;
    for i from 0 to n-1 do(
	partialSums := for i from 0 to m-1 list(sum(delete(0, flatten entries A_{i})));
	rowCheck#(#rowCheck) = (((unique sort partialSums) == {0,1}) or ((unique sort partialSums) == {0}) or ((unique sort partialSums) == {1}));
	);
    for i from 0 to m-1 do(
	partialSums := for i from 0 to n-1 list(sum(delete(0, flatten entries A^{i})));
	rowCheck#(#colCheck) = (((unique sort partialSums) == {0,1}) or ((unique sort partialSums) == {0}) or ((unique sort partialSums) == {1}));
	);
    (toList rowCheck == toList(#rowCheck:true)) and (toList colCheck == toList(#colCheck:true))
)

----------------------------------------
--Computes rank matrix of an ASM
--INPUT: an (n x n)- alternating sign matrix A OR a 1-line perm w
--OUTPUT: an (n x n) integer matrix of ranks of each entry
--Author: Yuyuan Luo
--TODO: add tests and documentation for this function
----------------------------------------

rankMatrix = method()
rankMatrix Matrix := Matrix => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
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
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    A := permToMatrix w;
    rankMatrix A
    )
-------------
--INPUT: a list w corresponding to a permutation in 1-line notation
    	--OR an alternating sign matrix A
--OUTPUT: a list of boxes in the Rothe diagram for A
--TODO: add documentation
-----------------------
rotheDiagram = method()
--this is a tidied version of code by Yuyuan Luo
rotheDiagram Matrix := List => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    n := numrows A;
    listEntries := flatten table(n,n, (i,j)->(i,j));
    ones := select(listEntries, i-> A_i == 1);
    seen := new MutableList;
    for one in ones do(
	for i from one_0 to n-1 do(
	    if (A_(i,one_1)==-1) then break;
	    seen#(#seen) = (i,one_1);
	    );
	for i from one_1 to n-1 do(
	    if A_(one_0,i)==-1 then break;
	    seen#(#seen) = (one_0,i);
	    );
	);
    seen = set unique toList seen;
    sort apply(toList((set listEntries) - seen), i-> (i_0+1,i_1+1))
    )

rotheDiagram List := List => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    A := permToMatrix w;
    rotheDiagram(A)
    )

-----------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
    	--OR an alternating sign matrix A
--OUTPUT: a list of essential boxes in the Rothe diagram for A
--TODO: add documentation
-----------------------
essentialBoxes = method()
essentialBoxes Matrix := List => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    boxes := rotheDiagram(A);
    badBoxes := apply(boxes, i->(positions(boxes,j->(j==(i_0,i_1+1)))|positions(boxes,j->(j==(i_0,i_1+1)))));
    essBoxes := positions(badBoxes, i-> i=={});
    boxes_essBoxes
    )
essentialBoxes List := List => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    essentialBoxes permToMatrix w
    )
----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: list of fulton generators for matrix determinantal ideal w
---------------------------------------
fultonGens = method()
fultonGens Matrix := List => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
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
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    A := permToMatrix w;
    fultonGens A
    )
----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: Schubert determinantal ideal for w
--WARNING: if you use the identity permutation your ring will be ZZ instead of Q and idk how to fix this
----------------------------------------
schubertDetIdeal = method()
schubertDetIdeal Matrix := Ideal => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    ideal fultonGens A
    )
schubertDetIdeal List := Ideal => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");    
    ideal fultonGens w
    );
----------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: single Grothendieck polynomials
--TODO: rename variables?
--WARNING: if you use the identity permutation your ring will be ZZ instead of Q and idk how to fix this
----------------------------
grothendieckPoly = method()
grothendieckPoly(Matrix):= (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    I := schubertDetIdeal A;
    Q := newRing(ring I, DegreeRank=> numcols A);
    numerator reduceHilbert hilbertSeries sub(I,Q)
    )
grothendieckPoly(List) := (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    I := schubertDetIdeal w;
    Q := newRing(ring I, DegreeRank=> #w);
    numerator reduceHilbert hilbertSeries sub(I,Q)
    )

--TODO: add tests
--TODO: double grothendieck. Problem: can't rename variables. use divided differences instead

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: ANTIdiagonal initial ideal of Schubert determinantal ideal for w
--WARNING: This method does not like the identity permutation
----------------------------------------
antiDiagInit = method()
antiDiagInit Matrix := monomialIdeal => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    monomialIdeal leadTerm schubertDetIdeal A
    );
antiDiagInit List := monomialIdeal => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    monomialIdeal leadTerm schubertDetIdeal w
    );
----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: diagonal initial ideal of Schubert determinantal ideal for w
--TODO: make diagRevlexInit function (potentially faster for tests)
--WARNING: This method does not like the identity permutation
----------------------------------------
diagLexInit = method()
diagLexInit Matrix := monomialIdeal => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    I:= schubertDetIdeal A;
    R:= newRing(ring I, MonomialOrder=>Lex); --making new ring with diagonal term order (lex suffices)
    monomialIdeal leadTerm sub(I,R)
    )
diagLexInit List := monomialIdeal => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
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
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    simplicialComplex antiDiagInit A
    );

subwordComplex List := simplicialComplex => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    simplicialComplex antiDiagInit w
    );

------------------------------------------
--INPUT: a nonempty list of equidimensional ASMs, presented as matrices
--OUTPUT: the minimal rank table, presented as a matrix
--TODO: tests, documentation
------------------------------------------

minimalRankTable = method()
minimalRankTable List := Matrix => (L) -> (
        if (#L == 0) then error("The input must be a nonempty list.");
        n := #(entries L#0);
        minimalRankMtx := mutableMatrix(ZZ,n,n);

        -- initialize the minimalRankMtx to something with big entries everywhere
        for i from 0 to n-1 do (
            for j from 0 to n-1 do (
                minimalRankMtx_(i,j) = n+1;
            );
        );

        -- comb through the list to get the minimal entrys
        for M in L do (
            listRankM := entries rankMatrix(M);
            if (#listRankM != n) then error ("The input must be a list of partial alternating sign matrices of the same size.");
            if not(isPartialASM(M)) then error("The input must be a list containing partial alternating sign matrices.");

            for i from 0 to n-1 do (
                for j from 0 to n-1 do (
                    minimalRankMtx_(i,j) = min {minimalRankMtx_(i,j), listRankM#i#j};
                );
            );
        );

        minimalRankMtx
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
      MatrixSchubert
    Headline
      a package for investigating matrix Schubert varieties and ASM varieties
    Description
      Text
       This package provides functions for constructing and investigating matrix Schubert varieties.
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
	  betti res antiDiagInit w	   
///

doc ///
    Key
	(isPartialASM, Matrix)
    	isPartialASM
    Headline
    	whether a matrix is a partial alternating sign matrix.
    Usage
    	isASM(M)
    Inputs
    	M:Matrix
    Description
    	Text
	 Given an integer matrix, checks that the matrix is a partial alternating sign matrix.
	 A partial alternating sign matrix is a matrix with entries in $\{-1,0,1\}$ such that:
	 
	     - The nonzero entries in each row and column alternate in sign,
	     
	     - Each row and column sums to $0$ or $1$, and
	     
	     - The first nonzero entry of any row or column is $1$.
	Example
    	    M = matrix{{0,0,1,0,0,0,0,0},{1,0,1,0,1,0,0,0},{0,0,0,1,-1,0,0,1},{0,0,1,-1,1,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,1,-1,1,0,0,0,0},{0,0,1,0,0,0,0,0}}
	    isPartialASM M
	    N = matrix{{0,-1,0,1,1},{1,-1,1,-1,1},{0,1,1,0,-1},{1,1,-1,1,-1},{-1,1,0,0,1}}
	    isPartialASM N

	    
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
  needsPackage "MatrixSchubert"
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
    matrix{{0,0,0,0,1,0,0,0},{0,0,1,0,-1,1,0,0},{0,0,0,1,0,0,0,0},{1,0,0,-1,1,-1,1,0},{0,1,-1,1,-1,1,0,0},{0,0,0,0,1,0,0,0},{0,0,1,0,0,0,0,0},{0,0,0,0,0,0,0,1}},
    matrix{{0,0,0},{0,1,0},{1,-1,0}},
    matrix{{1,0,0},{0,0,0}},
    matrix{{1,0,0},{0,0,1}},
    matrix{{0,1,0},{1,-1,0}},
    matrix{{0,1,0},{1,-1,0}},
    matrix{{0,0,1},{1,0,-1}},
    matrix{{0,0,1,0,0},{0,0,0,0,1},{0,0,0,0,0},{0,1,0,0,0}},
    matrix{{1,0,0,0},{0,0,1,0},{-1,1,0,0},{1,0,-1,1}}
    };
assert(apply(L,isPartialASM) == toList (#L:true))



T = {
    matrix{{-1}},
    matrix{{0,1,0},{1,0,1},{0,1,0}},
    matrix{{0,0,1,0,0,0,0,0},{1,0,1,0,1,0,0,0},{0,0,0,1,-1,0,0,1},{0,0,1,-1,1,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,1,-1,1,0,0,0,0},{0,0,1,0,0,0,0,0}}
    };
assert( apply(T, isPartialASM) == toList (#T:false))
///

TEST ///
--Example 2.1 in Weigandt "Prism Tableaux for ASMs"
A = matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}};
assert(isPartialASM A)
assert(sort essentialBoxes(A) == {(1,3),(2,1),(3,2)})

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
debug needsPackage "MatrixSchubert"

M = matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}}
w = {3,2,5,1,4}
permToMatrix w
isPartialASM M
rotheDiagram M
essentialBoxes M
grothendieckPoly M
netList fultonGens M
subwordComplex M
betti res diagLexInit M
betti res antiDiagInit M

------------------------------------
--Development Section
------------------------------------

restart
uninstallPackage "MatrixSchubert"
restart
installPackage "MatrixSchubert"
restart
needsPackage "MatrixSchubert"
elapsedTime check "MatrixSchubert"
viewHelp "MatrixSchubert"
