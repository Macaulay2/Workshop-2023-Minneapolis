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
    "minimalRankTable",
    "permLength",
    "augmentedRotheDiagram",
    "isPatternAvoiding",
    "isVexillary",
    "isRankTable"
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

------------------------------------
-- INPUT: A permutation in one-line notation
-- OUTPUT: The length of the permutation (number of inversions)
-- TODO: Add documentations + examples
------------------------------------

permLength = method()
permLength List := ZZ => w -> (
    if not (isPerm w) then error("the argument is not a permutation");
    k := 0;
    for i from 0 to #w - 2 do (
        for j from i + 1 to #w - 1 do (
            if w_i > w_j then k = k +1;
        ); 
    );
    k
)

--------------------------------
--checks if a permutation is pattern-avoiding
--INPUT: a permutation (in one-line notation), written as a list
--OUTPUT: whether the permutation avoid the pattern
--TODO: add documentation and tests
--------------------------------
isPatternAvoiding = method()
isPatternAvoiding (List,List) := Boolean => (perm, pattern) -> (
    -- Checks if the given permutation avoids the given pattern.
    -- Assume permutation is pattern-avoiding, break if not true
    isAvoiding := true;
    for idx in subsets(0..#perm-1, #pattern) do {
        sortedIdx := sort(idx);
        pairwiseComparison := apply(pattern_{0..#pattern-2}, pattern_{1..#pattern-1}, (i,j) -> perm#(sortedIdx#(i-1)) < perm#(sortedIdx#(j-1))); -- pairwise comparison of permutation according to pattern
        isAvoiding = not all(pairwiseComparison, i -> i == true); -- true if there was one inequality that failed, else all inequalities are true and so not pattern-avoiding
        if not isAvoiding then break;
    };
    isAvoiding
)

--------------------------------
--checks if a permutation is vexillary, i.e. 2143-avoiding
--INPUT: a permutation (one-line notation), written as a list
--OUTPUT: whether the permutation is vexillary
--TODO: add documentation and tests
--------------------------------
isVexillary = method()
isVexillary List := Boolean => (perm) -> (
    isPatternAvoiding(perm, {2,1,4,3})
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
--OUTPUT: A list of boxes in the Rothe diagram for A, with their corresponding rank values 
--TODO: add documentation + examples
-----------------------

augmentedRotheDiagram = method()
augmentedRotheDiagram List := List => w -> (
    L := rotheDiagram(w);
    R := rankMatrix(w);
    apply(L, (i,j) -> ((i,j), R_(i-1,j-1)))
)
augmentedRotheDiagram Matrix := List => w -> (
    L := rotheDiagram(w);
    R := rankMatrix(w);
    apply(L, (i,j) -> ((i,j), R_(i-1,j-1)))
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
--------------------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: Schubert determinantal ideal for w
--------------------------------------------
schubertDetIdeal = method()
schubertDetIdeal Matrix := Ideal => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    zMatrix := genMat(numrows A, numcols A); --generic matrix
    rankMat := rankMatrix A; --rank matrix for A
    essBoxes := essentialBoxes A;
    if essBoxes == {} then (
	R := ring zMatrix;
	return ideal(0_R)
	);
    zBoxes := apply(essBoxes, i-> flatten table(i_0,i_1, (j,k)->(j+1,k+1))); --smaller matrix indices for each essential box
    ranks := apply(essBoxes, i-> rankMat_(i_0-1,i_1-1)); --ranks for each essential box
    fultonGens := new MutableList;
    for box in essBoxes do (
    	pos := position(essBoxes, i-> i==box);
        fultonGens#(#fultonGens) = (minors(ranks_pos+1, zMatrix^{0..(box_0-1)}_{0..(box_1-1)}))_*;
        );
    return ideal(unique flatten toList fultonGens)
    );
schubertDetIdeal List := Ideal => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");    
    A := permToMatrix w;
    schubertDetIdeal A
    );

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: list of fulton generators for matrix determinantal ideal w
---------------------------------------
fultonGens = method()
fultonGens Matrix := List => (A) -> (
    (schubertDetIdeal A)_*
    );

fultonGens List := List => (w) -> (
    (schubertDetIdeal w)_*
    );

-*
----------------------------------------
--OLD CODE
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
--OLD CODE
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
*-

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

------------------------------------------
--INPUT: a square matrix M
--OUTPUT: whether M is a valid rank table.
--TODO: documentation, tests
------------------------------------------
isRankTable = method()
isRankTable Matrix := Boolean => (A) -> (
    AList := entries A;

    a := #AList;
    b := #(AList#0);
    if not(a == b) then return false;

    for i from 0 to a-1 do (    
        for j from 0 to a-1 do (
            if (i==0 and j==0 and not(AList#i#j == 0 or AList#i#j == 1)) then return false
            else if (i == 0 and j != 0 and (not(AList#i#j-AList#i#(j-1) == 0 or AList#i#j-AList#i#(j-1) == 1) or not(AList#i#j == 0 or AList#i#j == 1))) then return false
            else if (i != 0 and j == 0 and (not(AList#i#j-AList#(i-1)#j == 0 or AList#i#j-AList#(i-1)#j == 1) or not(AList#i#j == 0 or AList#i#j == 1))) then return false
            else if (i != 0 and j != 0 and (not((AList#i#j-AList#i#(j-1) == 0 or AList#i#j-AList#i#(j-1) == 1) or not(AList#i#j-AList#i#(j-1) == 0 or AList#i#j-AList#i#(j-1) == 1)))) then return false;
        );
    );

    true
);

-*
------------------------------------------
--INPUT: a rank table, presented as a matrix
--OUTPUT: an ASM corresponding to the rank table, presented as a matrix
------------------------------------------
rankTableToASM = method()
rankTableToASM Matrix := Matrix => (A) -> (
    n := #A;
    ASM :=  mutableMatrix(ZZ,n,n);

    -- find where all the ones are by checking whether it is larger than the entries above and to the left
    for i from 0 to n-1 do (
        for j from 0 to n-1 do (
            if (i == 0) then do (

            );
        );
    );
);*-


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
    	isPartialASM(M)
    Inputs
    	M:Matrix
    Description
    	Text
	 Given an integer matrix, checks that the matrix is a partial alternating sign matrix.
	 A partial alternating sign matrix is a matrix with entries in $\{-1,0,1\}$ such that:
	 
	     - The nonzero entries in each row and column alternate in sign,
	     
	     - Each row and column sums to $0$ or $1$, and
	     
	     - The first nonzero entry of any row or column (if there is one) is $1$.
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

doc ///
    Key
        (composePerms, List, List)
        composePerms
    Headline
        computes the composition of two permutations
    Usage
        composePerms(u,v)
    Inputs
        u:List
        v:List
    Description
        Text
            Computes the composition of two permutations, u and v, as u*v.
            Note that the permutations must be written as a list in one-line notation.
        Example
            u = {2,3,4,1}
            v = {4,3,2,1}
            composePerms(u,v)

            u = {1,2,3,4,5}
            v = {3,5,2,1,4}
            composePerms(u,v)

            u = {3,5,2,1,4}
            v = {1,2,3,4,5}
            composePerms(u,v)
///

doc ///
    Key
        (isPatternAvoiding, List, List)
        isPatternAvoiding
    Headline
        whether a permutation is pattern-avoiding, e.g. 2143-avoiding or 1432-avoiding
    Usage
        isPatternAvoiding(perm, pattern)
    Inputs
        perm:List
        pattern:List
    Description
        Text
            Given a permutation, checks if the permutation is pattern-avoiding, e.g. 2143-avoiding or 1432-avoiding.
            For example, a permutation w is 2143-avoiding if there does not exist indices i < j < k < l
            such that w_j < w_i < w_l < w_k.
        Example
            w = {7,2,5,8,1,3,6,4}
            pattern2143 = {2,1,4,3}
            isPatternAvoiding(w, pattern2143)

            v = {2,3,7,1,5,8,4,6}
            pattern1432 = {1,4,3,2}
            isPatternAvoiding(v, pattern1432)
///

doc ///
    Key
        (isVexillary, List)
        isVexillary
    Headline
        whether a permutation is vexillary, i.e. 2143-avoiding
    Usage
        isVexillary(perm)
    Inputs
        perm:List
    Description
        Text
            Given a permutation, checks if the permutation is vexillary, i.e. 2143-avoiding.
            A permutation w is 2143-avoiding if there do not exist indices i < j < k < l
            such that w_j < w_i < w_l < w_k.
        Example
            w = {7,2,5,8,1,3,6,4}
            isVexillary(w)

            v = {1,6,9,2,4,7,3,5,8}
            isVexillary(v)
///

-------------------------
-------------------------
--**TESTS SECTIONS**--
-------------------------
-------------------------

--detailed test for isPartialASM function
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

TEST ///
--composePerms
assert(composePerms({2,3,4,1}, {4,3,2,1}) == {1,4,3,2})
assert(composePerms({4,3,2,1}, {4,3,2,1}) == {1,2,3,4})
assert(composePerms({1,2,3,4,5}, {3,5,2,1,4}) == {3,5,2,1,4})
assert(composePerms({3,5,2,1,4}, {1,2,3,4,5}) == {3,5,2,1,4})
///

TEST ///
--isPatternAvoiding
assert(not isPatternAvoiding({2,3,7,1,5,8,4,6}, {1,4,3,2}));
assert(isPatternAvoiding({1,4,6,2,3,7,5}, {1,4,3,2}));

assert(not isPatternAvoiding({7,2,5,8,1,3,6,4}, {2,1,4,3}));
assert(isPatternAvoiding({1,6,9,2,4,7,3,5,8}, {2,1,4,3}));

--isVexillary
assert(not isVexillary({7,2,5,8,1,3,6,4}));
assert(isVexillary({1,6,9,2,4,7,3,5,8}));
///

end---------------------------------------------------------------


---------------------
--Ayah's sandbox
---------------------

restart
debug needsPackage "MatrixSchubert"

M = matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}}
fultonGens M
schubertDetIdeal M

w = {3,2,5,1,4}
schubertDetIdeal w

schubertDetIdeal {1,2,3}


permToMatrix w
isPartialASM M
rotheDiagram M
essentialBoxes M
grothendieckPoly M
netList fultonGens M
subwordComplex M
betti res diagLexInit M
betti res antiDiagInit M


A = matrix{{0,-1,0,1,1},{1,-1,1,-1,1},{0,1,1,0,-1},{1,1,-1,1,-1},{-1,1,0,0,1}}
n = numrows A
m = numcols A
rowCheck = new MutableList
colCheck = new MutableList
for i from 0 to n-1 do(
    partialSums = for i from 0 to m-1 list(sum(delete(0, flatten entries A_{i})));
    rowCheck#(#rowCheck) = (((unique sort partialSums) == {0,1}) or ((unique sort partialSums) == {0}) or ((unique sort partialSums) == {1}));
    );
for i from 0 to m-1 do(
    partialSums = for i from 0 to n-1 list(sum(delete(0, flatten entries A^{i})));
    rowCheck#(#colCheck) = (((unique sort partialSums) == {0,1}) or ((unique sort partialSums) == {0}) or ((unique sort partialSums) == {1}));
    );
(toList rowCheck == toList(#rowCheck:true)) and (toList colCheck == toList(#colCheck:true))
)
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

