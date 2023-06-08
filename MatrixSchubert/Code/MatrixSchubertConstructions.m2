
--------------------------------------------
--------------------------------------------
--**Constructing ASM Varieties**--
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
)


----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: ANTIdiagonal initial ideal of Schubert determinantal ideal for w
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


-*
--This will be the updated/faster antiDiagInit code
--but it's still under construction
antiDiagInit = method()
antiDiagInit Matrix := MonomialIdeal => (A) -> (
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
    antiDiagGens := new MutableList;
    for box in essBoxes do (
    	pos := position(essBoxes, i-> i==box);
	boxSubmatrix := zMatrix^{0..(box_0-1)}_{0..(box_1-1)};
	for x in subsets(numRows boxSubmatrix,) do (
	    for y in subsets(numCols boxSubmatrix do(
		    indicesList := 
        antiDiagGens#(#antiDiagGens) = apply((minors(ranks_pos+1, zMatrix^{0..(box_0-1)}_{0..(box_1-1)}))_*,);
        );
    return ideal(unique flatten toList antiDiagGens)
    );
antiDiagInit List := monomialIdeal => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    monomialIdeal leadTerm schubertDetIdeal w
    );
*-
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
    m := numcols A;
    listEntries := flatten table(n,m, (i,j)->(i,j));
    ones := select(listEntries, i-> A_i == 1);
    seen := new MutableList;
    for one in ones do(
	for i from one_0 to n-1 do( --death rays to the right
	    if (A_(i,one_1)==-1) then break;
	    seen#(#seen) = (i,one_1);
	    );
	for i from one_1 to m-1 do(
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

----------------------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation or an ASM ideal
--OUTPUT: whether or not R/I_A is CM
---------------------------------------
isSchubertCM = method()
isSchubertCM Matrix := Boolean => (A) -> (
    if not(isPartialASM A) then error("The input must be a partial alternating sign matrix or a permutation.");
    R:=ring(antiDiagInit A);
    codim(antiDiagInit A)==pdim(comodule (antiDiagInit A))
    );

isSchubertCM List := Boolean => (w) -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    print "We know from a theorem of Fulton that the quotient by any Schubert determinantal ideal is actually Cohen--Macaulay!";
    true
    );

isSchubertCM(matrix{{0,1,0},{1,-1,1},{0,1,0}})

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

        matrix minimalRankMtx
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

        matrix maximalRankMtx
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
    maxIdx := max((flatten entries vars ring I) / indexOfVariable // max);
    -- varWeights := (monoid ring I).Options.MonomialOrder#1#1;
    cycleDecomp := {};
    for primeComp in primeDecomp do {
        mons := sort(primeComp_*, mon -> ((indexOfVariable mon)_0+1)*maxIdx - (indexOfVariable mon)_1); --bad because variableIndex called twice; need decorated sort paradigm
        perms := apply(mons / indexOfVariable, perm -> toAntiDiagTrans(perm, maxIdx));
        cycleDecomp = append(cycleDecomp, fold(composePerms, perms));
    };
    unique cycleDecomp
)


-------------------------------------------
--INPUT: an ideal
--OUTPUT: whether the ideal is an intersection of Schubert determinantal ideals
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

-------------------------------------------
--INPUT: an ideal
--OUTPUT: whether the ideal is an ASM ideal
--TODO: docs and tests
--TODO: input validation/type checking
-------------------------------------------
isASMIdeal = method()
isASMIdeal Ideal := List => (I) -> (
    isASM := true;
    if (I == radical(I)) then {
        schubDecomp := schubertDecomposition I;
        if (isASM = I == intersect apply(schubDecomp/schubertDetIdeal, J -> sub(J, vars ring I))) then {
            permMatrices := (schubDecomp / permToMatrix);
            rankTable := rankTableFromMatrix matrix entrywiseMaxRankTable permMatrices;
            ASM := rankTableToASM matrix rankTable;
            ASMIdeal := schubertDetIdeal matrix ASM;
            isASM = I == sub(ASMIdeal, vars ring I);
            if isASM then I.cache.ASM = ASM;
        }
        else {
            isASM = false;
        }
    }
    else {
        isASM = false;
    };
    isASM
)

-------------------------------------------
--INPUT: an ideal
--OUTPUT: the ASM of an ideal
--WANRING: assumes the ideal is an ASM ideal
--TODO: docs and tests
--TODO: input validation/type checking
-------------------------------------------
getASM = method()
getASM Ideal := Matrix => (I) -> (
    I.cache.ASM
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
            else if (i != 0 and j != 0 and (not(AList#i#j-AList#i#(j-1) == 0 or AList#i#j-AList#i#(j-1) == 1) or not(AList#i#j-AList#(i-1)#j == 0 or AList#i#j-AList#(i-1)#j == 1))) then return false;
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

    matrix ASMret
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

    matrix rankTable
);

--------------------------------------------
-- INPUT: a list of permutations or ASMs
-- OUTPUT: the intersection of the ideals 
-- TODO: tests and documentation
--------------------------------------------
schubIntersect = method()
schubIntersect List := Ideal => (L) -> (
    if (#L == 0) then error("Please enter a nonempty list.");
    ll := L/schubertDetIdeal;
    intersect apply(ll, J -> sub(J, vars ring ll#0))
);

--------------------------------------------
-- INPUT: a list of permutations or ASMs
-- OUTPUT: the sum of the ideals 
-- TODO: tests and documentation
--------------------------------------------
schubAdd = method()
schubAdd List := Ideal => (L) -> (
    if (#L == 0) then error("Please enter a nonempty list.");
    listPermM := L / permToMatrix;
    rankM := entrywiseMinRankTable(listPermM);
    schubertDetIdeal rankTableToASM(rankM)
);