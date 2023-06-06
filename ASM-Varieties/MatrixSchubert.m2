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
	    "Posets",
            "MinimalPrimes"
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
    "schubertPoly",
    "doubleSchubertPoly",
    "entrywiseMinRankTable",
    "entrywiseMaxRankTable",
    "permLength",
    "augmentedRotheDiagram",
    "isPatternAvoiding",
    "isVexillary",
    "schubertDecomposition",
    "isIntersectionSchubIdeals",
    "rajCode",
    "rajIndex",
    "isMinRankTable",
    "Double",
    "rankTableToASM",
    "schubertReg",
    "rankTableFromMatrix",
    "schubertCodim"
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

isIdentity = method()
isIdentity List := Boolean => (w) -> (
    n := #w;
    w == toList(1..n)
    )

lastDescent = method()
lastDescent List := Boolean => (w) -> (
    if not (isPerm w) then error ("Expecting a permutation.");
    if isIdentity(w) then error ("Expecting a non-identity permutation.");
    n := #w;
   
    ans := -1;
    scan (reverse (0..n-2), i-> if w_i > w_(i+1) then (ans = i+1; break));
    ans
    )

permLength = method()
permLength List:= ZZ => (p) -> (
    if not(isPerm p) then error("The input must be a permutation.");
    l := 0;
    scan(#p, i->scan(i..#p-1, j ->(if p#i>p#j then l=l+1)));
     l)
 
swap = (L,i,j) -> (apply(#L, k-> if k!=i and k!=j then L_k
	                         else if k == i then L_j
				 else L_i))
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

--------------------------------
--auxiliary function for getting the index of a variable in a ring
--INPUT: an indexed variable
--OUTPUT: the index of the variable
--TODO: add docs and tests
-----------------------------------
variableIndex = method()
variableIndex IndexedVariable := Sequence => (elem) -> (
    --convert to string, parse, and select index, convert back
    value replace(".*_+", "", toString elem)
)
variableIndex RingElement := Sequence => (elem) -> (
    --convert to string, parse, and select index, convert back
    value replace(".*_+", "", toString elem)
)

-----------------------------------------------
-----------------------------------------------
--**Useful functions for permutations**
-----------------------------------------------
-----------------------------------------------


------------------------------------
--INPUT: A transposition in cycle notation, and the n for which to regard perm 
--       as an element of S_n
--OUTPUT: the transposition in one-line notation
--TODO: docs and tests
------------------------------------
toOneLineNotation = method()
toOneLineNotation (Sequence, ZZ) := List => (perm, maxIdx) -> (
    switch(perm_0-1, perm_1-1, toList(1..maxIdx))
)
toOneLineNotation (List, ZZ) := List => (perm, maxIdx) -> (
    switch(perm_0-1, perm_1-1, toList(1..maxIdx))
)

------------------------------------
--INPUT: An index (i,j)
--OUTPUT: the corresponding transposition according to antidiagonal term order
--TODO: docs and tests
------------------------------------
toAntiDiagTrans = method()
toAntiDiagTrans (Sequence, ZZ) := List => (idx, maxIdx) -> (
    transposition := (sum(toList idx)-1, sum(toList idx));
    toOneLineNotation(transposition, maxIdx)
)
toAntiDiagTrans (List, ZZ) := List => (idx, maxIdx) -> (
    transposition := {sum(toList idx)-1, sum(toList idx)};
    toOneLineNotation(transposition, maxIdx)
)

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

-- ------------------------------------
-- -- INPUT: A permutation in one-line notation
-- -- OUTPUT: The length of the permutation (number of inversions)
-- -- TODO: Add documentations + examples
-- ------------------------------------

-- permLength = method()
-- permLength List := ZZ => w -> (
--     if not (isPerm w) then error("the argument is not a permutation");
--     k := 0;
--     for i from 0 to #w - 2 do (
--         for j from i + 1 to #w - 1 do (
--             if w_i > w_j then k = k +1;
--         ); 
--     );
--     k
-- )

--------------------------------
--checks if a permutation is pattern-avoiding
--INPUT: a permutation (in one-line notation), written as a list
--OUTPUT: whether the permutation avoid the pattern
--TODO: input validation/type checking
--------------------------------
isPatternAvoiding = method()
isPatternAvoiding (List,List) := Boolean => (perm, pattern) -> (
    --input validation
    if not (isPerm perm) then error(toString perm | " is not a permutation.");
    --assume permutation is pattern-avoiding, break if not true
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
--TODO: input validation/type checking
--------------------------------
isVexillary = method()
isVexillary List := Boolean => (perm) -> (
    --input validation
    if not (isPerm perm) then error(toString perm | " is not a permutation.");
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
    n:=numrows(A);
    m:=numcols(A);
for i from 0 to n-1 do (
    rowPartialSum := accumulate(plus,{0}|(flatten entries(A^{i})));
    if (not(isSubset(sort unique rowPartialSum,{0,1}))) then return false;
    );
for i from 0 to m-1 do (
    colPartialSum := accumulate(plus,{0}|(flatten entries(A_{i})));
    if (not(isSubset(sort unique colPartialSum, {0,1}))) then return false;
    );
return true
); 

----------------------------------------
--Computes rank matrix of an ASM
--INPUT: an (n x n)- alternating sign matrix A OR a 1-line perm w
--OUTPUT: an (n x n) integer matrix of ranks of each entry
--Author: Yuyuan Luo
--TODO: add tests for this function
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
--TODO:
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



schubertPolyHelper = method(Options=>{Double=>false})
schubertPolyHelper(List, Ring) := opts->(w, Q) -> (
    isDouble := opts.Double;
    n := #w;
    if not (isPerm w) then error ("The input must be a permutation matrix.");
   -- Q := QQ[x_1..x_n];
    if (isIdentity w) then return 1;
    r := lastDescent(w) - 1;
    --print(r);
    --print(reverse(r+1..n-1));
    s := -1;
    scan(reverse(r+1..n-1), i-> if w_i < w_r then (s = i; break));
    --print({r,s});
    v := swap(w,r,s);
    previnds := select(0..r-1, q-> permLength(swap(v,q,r))==permLength(v)+1);
    us := apply(previnds, i-> swap(v,i,r));
    if not isDouble then
        sum(toList(apply(us, u->schubertPolyHelper(u,Q,Double=>isDouble))))+ Q_r * schubertPolyHelper(v,Q,Double=>isDouble)
    else 
        sum(toList(apply(us, u->schubertPolyHelper(u,Q,Double=>isDouble))))+ (Q_r-Q_(n-1+v_r)) * schubertPolyHelper(v,Q,Double=>isDouble)
    )

schubertPoly = method()
schubertPoly(List) := (w) -> (
    n := #w;
    x := local x;
    Q := QQ[x_1..x_n];
    schubertPolyHelper(w, Q, Double=>false)
    )

doubleSchubertPoly = method()
doubleSchubertPoly(List) := (w) -> (
    n := #w;
    x := local x;
    y := local y;
    Q := QQ[x_1..x_n,y_1..y_n];
    schubertPolyHelper(w, Q, Double=>true)
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
--INPUT: a list w corresponding to a permutation in 1-line notation or a partial ASM
--OUTPUT: the Castlenuovo-Mumford regularity of I_A or I_w
----------------------------------------
schubertReg = method()
schubertReg Matrix := ZZ => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    regularity(antiDiagInit A)
    );
schubertReg List := ZZ => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
     regularity(antiDiagInit w)
    );

schubertCodim = method() 
schubertCodim Matrix := ZZ => A -> (
    if not (isPartialASM A) then error("The input must be a partial alternating sign matrix");
    codim antiDiagInit A
)
schubertCodim List := ZZ => (w) -> (
    if not (isPerm w) then error("The input must be a permutation in one line notation");
    permLength w
)

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

subwordComplex List := simplicialComplex => (w) -> (
    if not(isPerm w) then error("The input must be a permutation.");
    simplicialComplex antiDiagInit w
    );

------------------------------------------
--INPUT: a nonempty list of equidimensional ASMs, presented as matrices
--OUTPUT: the minimal rank table, presented as a matrix
--TODO: tests, documentation
------------------------------------------

entrywiseMinRankTable = method()
entrywiseMinRankTable List := Matrix => (L) -> (
        if (#L == 0) then error("The input must be a nonempty list.");
        n := #(entries L#0);
        minimalRankMtx := mutableMatrix(ZZ,n,n);

        -- initialize the minimalRankMtx to something with big entries everywhere
        for i from 0 to n-1 do (
            for j from 0 to n-1 do (
                minimalRankMtx_(i,j) = n+1;
            );
        );

        -- comb through the list to get the minimal entries
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
--INPUT: a nonempty list of equidimensional ASMs, presented as matrices
--OUTPUT: the maximal rank table, presented as a matrix
--TODO: tests, documentation
------------------------------------------

entrywiseMaxRankTable = method()
entrywiseMaxRankTable List := Matrix => (L) -> (
        if (#L == 0) then error("The input must be a nonempty list.");
        n := #(entries L#0);
        maximalRankMtx := mutableMatrix(ZZ,n,n);


        -- comb through the list to get the maximal entries
        for M in L do (
            listRankM := entries rankMatrix(M);
            if (#listRankM != n) then error ("The input must be a list of partial alternating sign matrices of the same size.");
            if not(isPartialASM(M)) then error("The input must be a list containing partial alternating sign matrices.");

            for i from 0 to n-1 do (
                for j from 0 to n-1 do (
                    maximalRankMtx_(i,j) = max {maximalRankMtx_(i,j), listRankM#i#j};
                );
            );
        );

        maximalRankMtx
    );

-------------------------------------------
--INPUT: an ASM ideal
--OUTPUT: the primary decomposition of the ASM ideal
--TODO: docs and tests
--TODO: input validation/type checking
-------------------------------------------
schubertDecomposition = method()
schubertDecomposition Ideal := List => (I) -> (
    primeDecomp := decompose ideal leadTerm I;
    maxIdx := max((flatten entries vars ring I) / variableIndex // max);
    -- varWeights := (monoid ring I).Options.MonomialOrder#1#1;
    cycleDecomp := {};
    for primeComp in primeDecomp do {
        mons := sort(primeComp_*, mon -> ((variableIndex mon)_0+1)*maxIdx - (variableIndex mon)_1); --bad because variableIndex called twice; need decorated sort paradigm
        perms := apply(mons / variableIndex, perm -> toAntiDiagTrans(perm, maxIdx));
        cycleDecomp = append(cycleDecomp, fold(composePerms, perms));
    };
    unique cycleDecomp
)


-------------------------------------------
--INPUT: an ideal
--OUTPUT: whether the ideal is an ASM ideal
--WARNING: Might not be right depending on if {2,1,6,3,5,4} is ASM.
--         If it is, then buggy; else, we are good.
--TODO: docs and tests
--TODO: input validation/type checking
-------------------------------------------
isIntersectionSchubIdeals = method()
isIntersectionSchubIdeals Ideal := List => (I) -> (
    isIntersection := true;
    if (I == radical(I)) then {
        schubDecomp := schubertDecomposition I;
        isIntersection = I == intersect apply(schubDecomp/schubertDetIdeal, J -> sub(J, vars ring I));
    }
    else {
        isIntersection = false;
    };
    isIntersection
)
------------------------------------------
--INPUT: a square matrix M
--OUTPUT: whether M is a valid rank table.
--TODO: documentation, tests
------------------------------------------
isMinRankTable = method()
isMinRankTable Matrix := Boolean => (A) -> (
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


------------------------------------------
--INPUT: a rank table, presented as a matrix
--OUTPUT: an ASM corresponding to the rank table, presented as a matrix
--TODO: documentation and tests
------------------------------------------
rankTableToASM = method()
rankTableToASM Matrix := Matrix => (A) -> (
    if not(isMinRankTable(A)) then error("The inputted matrix is not a valid minimal rank table.");
    AList := entries A;
    n := #AList;
    ASMret :=  mutableMatrix(ZZ,n,n);
    for i from 0 to n-1 do (
        for j from 0 to n-1 do (
            if (i == 0 and j == 0) then (
                if (AList#0#0 == 1) then (A_(0,0) = 1;);
            )
            else if (i == 0) then (
                if (AList#i#j == 1 and AList#i#(j-1)==0) then (ASMret_(i,j) = 1;);
            )
            else if (j == 0) then (
                if (AList#i#j == 1 and AList#(i-1)#j==0) then (ASMret_(i,j) = 1;);
            )
            else (
                if (AList#i#j - AList#i#(j-1) == 1 and AList#i#j - AList#(i-1)#j == 1 and AList#(i-1)#j == AList#(i-1)#(j-1)) then (ASMret_(i,j) = 1;)
                else if (AList#i#j == AList#i#(j-1) and AList#i#j == AList#(i-1)#j and AList#i#j > AList#(i-1)#(j-1)) then (ASMret_(i,j) = -1;);
            );
        );
    );

    ASMret
);

--------------------------------------------
-- INPUT: an integer matrix M where the entries are at least 0
-- OUTPUT: the minimal rank table associated to M representing an ASM 
-- TODO: tests and documentation
--------------------------------------------
rankTableFromMatrix = method()
rankTableFromMatrix Matrix := Matrix => (A) -> (
    AList := entries A;
    n := #AList;
    rankTable := mutableMatrix(ZZ,n,n);
    if not(#(AList#0) == n) then error("Must be a square matrix.");

    for i from 0 to n-1 do (
        for j from 0 to n-1 do(
            if not(ring(AList#i#j) === ZZ) then error("Must be an integer matrix.");
            if (AList#i#j < 0) then error("Must be a matrix with nonnegative entries.");
            if (i == 0 and j == 0) then (
                rankTable_(n-1,n-1) = min(n,AList#(n-1)#(n-1));
            )
            else if (i == 0) then (
                rankTable_(n-1-i,n-1-j) = min(n-i,n-j,AList#(n-1-i)#(n-1-j),rankTable_(n-1-i,n-j));
            )
            else if (j == 0) then (
                rankTable_(n-1-i,n-1-j) = min(n-i,n-j,AList#(n-1-i)#(n-1-j),rankTable_(n-i,n-j-1));
            )
            else (
                rankTable_(n-1-i,n-1-j) = min(n-i,n-j,AList#(n-1-i)#(n-1-j),rankTable_(n-1-i,n-j),rankTable_(n-i,n-j-1));
            );
        );
    );

    for i from 0 to n-1 do (
        for j from 0 to n-1 do(
            if (i == 0 and j == 0) then (
                rankTable_(i,j) = min(1,rankTable_(i,j));
            )
            else if (i == 0) then (
                rankTable_(i,j) = min(rankTable_(i,j), rankTable_(i,j-1)+1);
            )
            else if (j == 0) then (
                rankTable_(i,j) = min(rankTable_(i,j), rankTable_(i-1,j)+1);
            )
            else (
                rankTable_(i,j) = min(rankTable_(i,j), rankTable_(i,j-1)+1, rankTable_(i-1,j)+1);
            );
        );
    );

    rankTable
);


----------------------------------------
-- Part 2. Invariants of ASM Varieties
----------------------------------------

------------------------------------------
--INPUT: lengthIncrSeq, takes a permutation in one line notation
--OUTPUT: returns the length of the longest consecutive permutation
--        which starts at the beginning of the permutation
--TO DO: Thoroughly test and document
------------------------------------------

lengthIncrSubset = (w) -> (
   
   if (w == {}) then return 0;
   
   preVal := w_0;
   for i from 1 to #w-1 do (
       if (preVal > w_i) then return i;
       preVal = w_i;
   );
   return #w;
);


------------------------------------------
--INPUT: rajCode, takes a permutation in one line notation
--OUTPUT: returns the rajCode of the permutation
------------------------------------------

rajCode = method()
rajCode List := ZZ => (w) -> (

    if not (isPerm w) then error ("Expecting a permutation.");
   
    rajCodeVec := {};
    for k from 0 to #w-2 do (
	maxLengthIncr := 1;
	fVal := w_k;
	subPerm := w_{k+1..#w-1};
	
	for l in delete({},subsets(subPerm)) do (
	    testPerm := {fVal} | l;
	    maxLengthIncr = max(maxLengthIncr,lengthIncrSubset(testPerm));
	);
    	
	rajCodeVec := rajCodeVec | {#subPerm+1 - maxLengthIncr};
    );
    return rajCodeVec;
);


------------------------------------------
--INPUT: rajIndex, takes a permutation in one line notation
--OUTPUT: returns the rajIndex of the permutation
------------------------------------------
rajIndex = method()
rajIndex List := ZZ => (w) -> (
  
    if not (isPerm w) then error ("Expecting a permutation.");
    
    return sum rajCode w;
    
);



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
	(permLength, List)
        permLength
    Headline
    	to find the length of a permutation in 1-line notation.
    Usage
        permLength(w)
    Inputs
    	w:List
    Description
    	Text
	 Given a permutation in 1-line notation returns the Coxeter length of the permutation.
	Example
    	    w = {2,5,4,1,3}
	    permLength(w)

	    
///

doc ///
    Key
        (rotheDiagram, List)
	(rotheDiagram, Matrix)
    	rotheDiagram
    Headline
    	to find the Rothe diagram of a partial alternating sign matrix
    Usage
    	rotheDiagram(w)
	rotheDiagram(M)
    Inputs
    	w:List
	    or {\tt M} is a @TO Matrix@
    Description
    	Text
	 Given a permutation in 1-line notation or a partial alternating sign matrix returns the Rothe diagram.
	Example
    	    w = {2,5,4,1,3}
	    rotheDiagram(w)
	    M = matrix{{0,1,0},{1,-1,0},{0,0,0}}
	    rotheDiagram(M)

	    
///

doc ///
    Key
        (augmentedRotheDiagram, List)
	(augmentedRotheDiagram, Matrix)
    	augmentedRotheDiagram
    Headline
    	to find the Rothe diagram of a partial alternating sign matrix together with the rank conditions determining the alternating sign matrix variety
    Usage
    	augmentedRotheDiagram(w)
	augmentedRotheDiagram(M)
    Inputs
    	w:List
	    or {\tt M} is a @TO Matrix@
    Description
    	Text
	 Given a permutation in 1-line notation or a partial alternating sign matrix returns list of entries of Rothe diagram with the ranks of each entry.
	Example
    	    w = {2,5,4,1,3}
	    augmentedRotheDiagram(w)
	    M = matrix{{0,1,0},{1,-1,0},{0,0,0}}
	    augmentedRotheDiagram(M)

	    
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
    	(schubertDetIdeal, Matrix)
    	schubertDetIdeal
    Headline
    	Computes Schubert determinantal ideal for a given permutation.
    Usage
    	schubertDetIdeal(M)
    Inputs
    	M:Matrix
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

doc ///
    Key
        (schubertDecomposition, Ideal)
        schubertDecomposition
    Headline
        finds the decomposition of an ASM ideal into Schubert determinantal ideals
    Usage
        schubertDecomposition(I)
    Inputs
        I:Ideal
    Description
        Text
            Given an ASM ideal, it can be decomposed into Schubert determinantal ideals
            as I = I_{w_1} \cap ... \cap I_{w_k}, where the w_i are permutations.
            As output, each element in the list is the permutation associated 
            to a prime component in the Schubert decomposition of the antidiagonal 
            initial ideal of I.
///

doc ///
    Key
        (isIntersectionSchubIdeals, Ideal)
        isIntersectionSchubIdeals
    Headline
        whether an ideal is ASM
    Usage
        isIntersectionSchubIdeals(I)
    Inputs
        I:Ideal
    Description
        Text
            An ideal I is ASM if I is radical and I = I_{w_1} \cap ... \cap I_{w_k},
            where the I_{w_i} are Schubert determinantal ideals.
///

doc ///
    Key
    	rankMatrix
	(rankMatrix, Matrix)
	(rankMatrix, List)
    Headline
    	Computes a matrix of rank conditions that determines a Schubert determinantal ideal or, more generally, an alternating sign matrix ideal.
    Usage
    	rankMatrix(M)
	rankMatrix(w)
    Inputs
    	M:Matrix
	w:List
    Description
    	Text
	 Given an alternating sign matrix or a permutation in 1-line notation, 
	 outputs the matrix of rank condition associated to that alternating sign matrix or permutation.
	Example
	 rankMatrix({1,3,2})
	 rankMatrix(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

///

doc ///
    Key
    	essentialBoxes
	(essentialBoxes, Matrix)
	(essentialBoxes, List)
    Headline
    	Computes a list of the essential boxes in the Rothe Diagram for an alternating sign matrix M or a permutation in 1-line notation.
    Usage
    	essentialBoxes(M)
	essentialBoxes(w)
    Inputs
    	M:Matrix
	w:List
    Description
    	Text
	 Given an alternating sign matrix or a permutation in 1-line notation, 
	 outputs a list of the essential boxes in the Rothe diagram for that matrix. 
	Example
	 essentialBoxes({1,3,2})
	 essentialBoxes(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

///

doc ///
    Key
        schubertCodim
        (schubertCodim, Matrix)
        (schubertCodim, List)
    Headline
        Computes the codimension of a schubert determinantal ideal
    Usage
        schubertCodim(M)
        schubertCodim(w)
    Inputs
        w:List
	        or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, outputs the codimension of the corresponding schubert determinantal ideal.
        Example
            schubertCodim(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
            schubertCodim({1,3,2})
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
    matrix{{0,0,1,0,0},{0,0,0,0,1},{0,0,0,0,0},{0,1,0,0,0}}
    };
assert(apply(L,isPartialASM) == toList (#L:true))



T = {
    matrix{{-1}},
    matrix{{0,1,0},{1,0,1},{0,1,0}},
    matrix{{0,0,1,0,0,0,0,0},{1,0,1,0,1,0,0,0},{0,0,0,1,-1,0,0,1},{0,0,1,-1,1,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,1,-1,1,0,0,0,0},{0,0,1,0,0,0,0,0}},
    matrix{{1,0,0,0},{0,0,1,0},{-1,1,0,0},{1,0,-1,1}}
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

TEST /// 
-- permLength 

assert(permLength {1} == 0)
assert(permLength {1,2} == 0)
assert(permLength {3,2,1} == 3)
assert(permLength {2,1,3} == 1)
assert(permLength {8,7,6,5,4,3,2,1} == 28)
///

TEST ///
-- augmentedRotheDiagram 

assert(sort augmentedRotheDiagram {2,1,5,4,3} == sort {((1,1),0), ((3,3),2),((3,4),2), ((4,3),2)})
assert(sort augmentedRotheDiagram matrix{{0,1,0},{1,-1,1},{0,1,0}} == sort{((1,1),0), ((2,2),1)})
assert (sort augmentedRotheDiagram matrix {{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}} == sort {((1,1),0),((1,2),0),((4,3),2),((3,3),2)})
///


TEST ///
-- essentialBoxes 

assert(essentialBoxes({2,1,6,3,5,4 })== {(1, 1), (3, 5), (5, 4)})
assert(essentialBoxes matrix {{0,1,0,0,0,0},{1,0,0,0,0,0},{0,0,0,0,0,1},{0,0,1,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0}} == {(1, 1), (3, 5), (5, 4)})
///


TEST ///
-- schubertCodim
L = {
    {1},
    {2,1},
    matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}},
    matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}}
}
assert all (L, i -> schubertCodim i == codim schubertDetIdeal i)
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
