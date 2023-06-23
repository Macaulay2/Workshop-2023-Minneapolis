-----------------------------------------------
-----------------------------------------------
--**Useful functions for permutations**
-----------------------------------------------
-----------------------------------------------

-----------------------------------------------------------
--tests if a list is a permutation in 1-line notation
--INPUT: A list of integers w
--OUTPUT: TRUE if w is a permutation; else FALSE
-----------------------------------------------------------
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


-----------------------------------------------------------
--tests of a permutation is the identity permutation
--INPUT: a list of integers w
--OUTPUT: TRUE if w is the identity perm; else FALSE
-----------------------------------------------------------

isIdentity = method()
isIdentity List := Boolean => (w) -> (
    n := #w;
    w == toList(1..n)
)

-----------------------------------------------------------
--finds index of last descent in a permutation
--INPUT: a permutation in 1-line notation
--OUPTUT: last index where w_i > w_(i+1)
-----------------------------------------------------------


lastDescent = method()
lastDescent List := ZZ => (w) -> (
    if not (isPerm w) then error ("Expecting a permutation.");
    if isIdentity(w) then error ("Expecting a non-identity permutation.");
    n := #w;
   
    ans := -1;
    scan (reverse (0..n-2), i -> if w_i > w_(i+1) then (ans = i+1; break));
    ans
)

-----------------------------------------------------------
--find index of first descent in a permutation
--INPUT: a permutation in 1-line notation
--OUTPUT: first index where w_i > w_(i+1)
-----------------------------------------------------------


firstDescent = method()
firstDescent List := ZZ => (w) -> (
    if not (isPerm w) then error ("Expecting a permutation.");
    if isIdentity(w) then error ("Expecting a non-identity permutation.");
    n := #w;
   
    ans := -1;
    scan ((0..n-2), i-> if w_i > w_(i+1) then (ans = i+1; break));
    ans
)

-----------------------------------------------------------
--compute the length of a permutation
-----------------------------------------------------------

permLength = method()
permLength List:= ZZ => (p) -> (
    if not(isPerm p) then error("The input must be a permutation.");
    l := 0;
    scan(#p, i -> scan(i..#p-1, j -> (if p#i > p#j then l = l+1)));
    l
)
 
-----------------------------------------------------------
--creates simple transposition in a permtuation
-----------------------------------------------------------

swap = (L,i,j) -> (
    apply(#L, k -> if k != i and k != j then L_k
	               else if k == i then L_j
				   else L_i)
)

-----------------------------------------------------------
--computes inverse permutation
-----------------------------------------------------------

inverseOf = method()
inverseOf List := List => (w) -> (
    if not (isPerm w) then error ("Expecting a permutation.");
    winv := new MutableList;
    scan(#w, i-> winv#(w_i-1)=i+1);
    toList winv
    )

-----------------------------------------------------------
--creates longest word/longest permutation {n,n-1,...2,1}
-----------------------------------------------------------

longestPerm = method()
longestPerm ZZ := List => (n) -> (
    if n < 1 then error ("Expecting a positive integer.");
    toList (reverse (1..n))
    )

-----------------------------------------------------------
--INPUT: a permutation as a list of integers w
--OUPT: list of simple transpositions giving a reduced word for w
-----------------------------------------------------------

getOneReducedWord = method()
getOneReducedWord List := List => (w) -> (
    if not (isPerm w) then error ("Expecting a permutation.");
    if isIdentity(w) then  {}
    else (
        i := firstDescent(w);
        newPerm := swap(w, i-1, i);
        getOneReducedWord(newPerm)) | {i}
    )


------------------------------------
--INPUT: A transposition in cycle notation, and the n for which to regard perm 
--       as an element of S_n
--OUTPUT: the transposition in one-line notation
--TODO: docs and tests
------------------------------------
toOneLineNotation = method()
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
    transposition := {sum(toList idx)-1, sum(toList idx)};
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
    if not (isPerm v) then error("The second argument is not a permutation.");
    if not (#v==#u) then error("Expected permutations of the same length.");
    u0 := apply(u, i->i-1);
    v0 := apply(v, i->i-1);
    apply(u0_v0, i-> i+1)
)

--------------------------------
--checks if a permutation is pattern-avoiding
--INPUT: a permutation (in 1-line notation), written as a list
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
--INPUT: a permutation (1-line notation), written as a list
--OUTPUT: whether the permutation is vexillary
--TODO: input validation/type checking
--------------------------------
isVexillary = method()
isVexillary List := Boolean => (perm) -> (
    if not (isPerm perm) then error(toString perm | " is not a permutation.");
    isPatternAvoiding(perm, {2,1,4,3})
)

--------------------------------
--checks if a permutation avoids all of the given patterns, i.e. 2143-avoiding
--INPUT: a permutation (1-line notation), written as a list, and a list of lists (patterns)
--OUTPUT: whether the permutation avoids all of the patterns
--TODO: input validation/type checking
--------------------------------
avoidsAllPatterns = method()
avoidsAllPatterns (List, List) := Boolean => (perm, patterns) -> (
    all(patterns, pattern -> isPatternAvoiding(perm, pattern))
)

--------------------------------
--checks if a permutation is Cartwright-Sturmfels
--INPUT: a permutation (1-line notation)
--OUTPUT: whether the permutation avoids all of the patterns
--TODO: input validation/type checking
--------------------------------
isCartwrightSturmfels = method()
isCartwrightSturmfels List := Boolean => (perm) -> (
    patterns := {{1,2,5,4,3}, 
                 {1,3,2,5,4}, 
                 {1,3,5,2,4}, 
                 {1,3,5,4,2}, 
                 {2,1,5,4,3}, 
                 {1,2,5,3,6,4}, 
                 {1,2,5,6,3,4}, 
                 {2,1,5,3,6,4},
                 {2,1,5,6,3,4},
                 {3,1,5,2,6,4},
                 {3,1,5,6,2,4},
                 {3,1,5,6,4,2}};
    all(patterns, pattern -> isPatternAvoiding(perm, pattern))
)

--------------------------------
--checks if a permutation is CDG
--INPUT: a permutation (1-line notation)
--OUTPUT: whether the permutation avoids all of the patterns
--TODO: input validation/type checking
--------------------------------
isCDG = method()
isCDG List := Boolean => (perm) -> (
    patterns := {{1,3,2,5,4},
                 {2,1,5,4,3},
                 {2,1,4,6,3,5},
                 {2,1,5,3,6,4},
                 {2,1,5,6,3,4},
                 {2,4,1,6,3,5},
                 {3,1,5,2,6,4},
                 {4,2,6,1,7,3,5}};
    all(patterns, pattern -> isPatternAvoiding(perm, pattern))
)

------------------------------------------
--INPUT: longestIncrSeq, takes a previous value, a previous size, and a permutation in one line notation
--OUTPUT: returns the length of the longest consecutive permutation plus the previous size
--        that starts at the beginning of the permutation and has the elements of the sequence larger than the previous value
--TODO: Document
------------------------------------------
longestIncrSeq = method()
longestIncrSeq (ZZ,ZZ,List) := List => memoize ((preVal,prevSZ,w) -> (
    if w == {} then return prevSZ;
    
    longestSZ := prevSZ;
    currSZ := prevSZ;
    currVal := preVal;
    for i from 0 to #w-1 do (
    	currVal = w_i;
	if (currVal < preVal) then currSZ = longestIncrSeq(preVal, prevSZ, w_{i+1..#w-1})
	else currSZ = longestIncrSeq(currVal, prevSZ+1, w_{i+1..#w-1});
	longestSZ = max(longestSZ, currSZ);
    );
    return longestSZ;
))

------------------------------------------
--INPUT: rajCode, takes a permutation in one line notation
--OUTPUT: returns the rajCode of the permutation
------------------------------------------
rajCode = method()
rajCode List := ZZ => (w) -> (
    if not (isPerm w) then error ("Expecting a permutation.");
    rajCodeVec := {};
    for k from 0 to #w-1 do (
    	rajCodeVec = append(rajCodeVec, (#w-k) - longestIncrSeq(w_k, 1, w_{k+1..#w-1}));
    );
    return rajCodeVec;
)

------------------------------------------
--INPUT: rajIndex, takes a permutation in one line notation
--OUTPUT: returns the rajIndex of the permutation
--TO DO: document
------------------------------------------
rajIndex = method()
rajIndex List := ZZ => (w) -> ( 
    if not (isPerm w) then error ("Expecting a permutation.");
    return sum rajCode w;
)

----------------------------
--INPUT: a list w corresponding to a permutation in 1-line notation
--OUTPUT: single Grothendieck polynomials
--TODO: rename variables?
--This function is now correct, but it's slow
--potentially we could rename it to "Kpolynomial" or something
    --then we could allow it to take ASMs as input
----------------------------
grothendieckPoly = method(Options=>{Algorithm=>"DividedDifference"})
grothendieckPoly(List) := opts -> w -> (
    if not(isPerm w) then error("The input must be a partial alternating sign matrix or a permutation.");
    if opts.Algorithm == "Degree" then (
        I := schubertDetIdeal w;
        R := ring I;
        kk := coefficientRing R;
        possibleDegs := apply(#w, i-> toList insert(i,1,(#w-1):0));
        degs := splice apply(possibleDegs, i->#w:i);
        Q := kk[R_*, Degrees => degs];
        numerator hilbertSeries sub(I,Q)
    )
    else if opts.Algorithm == "DividedDifference" then (
        n := #w;
        x := local x;
        Q = QQ[x_1..x_n];
        polyByDividedDifference(w, Q, PolyType => "Grothendieck")        	
    )
    else error("Invalid option for Algorithm.")
)


schubertPolyHelper = method(Options=>{Double=>false})
schubertPolyHelper(List, PolynomialRing) := opts -> (w, Q) -> (
    isDouble := opts.Double;
    n := #w;
    if not (isPerm w) then error ("The input must be a permutation matrix.");
    if (isIdentity w) then return 1;
    r := lastDescent(w) - 1;
    s := -1;
    scan(reverse(r+1..n-1), i -> if w_i < w_r then (s = i; break));
    v := swap(w,r,s);
    previnds := select(0..r-1, q -> permLength(swap(v,q,r)) == permLength(v)+1);
    us := apply(previnds, i -> swap(v, i, r));
    if not isDouble then
        sum(toList(apply(us, u -> schubertPolyHelper(u, Q, Double=>isDouble)))) + Q_r * schubertPolyHelper(v, Q, Double=>isDouble)
    else 
        sum(toList(apply(us, u -> schubertPolyHelper(u, Q, Double=>isDouble)))) + (Q_r - Q_(n-1+v_r)) * schubertPolyHelper(v, Q, Double=>isDouble)
)

schubertPoly = method(Options=>{Algorithm=>"DividedDifference"})
schubertPoly(List) := opts-> (w) -> (
    
    n := #w;
    x := local x;
    Q := QQ[x_1..x_n];
    if opts.Algorithm == "Transition" then 
        schubertPolyHelper(w, Q, Double=>false)
    else if opts.Algorithm == "DividedDifference" then
        polyByDividedDifference(w,Q)
    else error("Invalid option for Algorithm.")
    
)

doubleSchubertPoly = method()
doubleSchubertPoly(List) := (w) -> (
    n := #w;
    x := local x;
    y := local y;
    Q := QQ[x_1..x_n,y_1..y_n];
    schubertPolyHelper(w, Q, Double=>true)
)



dividedDifference = method(Options=>{Operator=>null})
dividedDifference (RingElement, ZZ) := opts-> (f,i) -> (
    Q:= ring f;
    sf := sub(f, {Q_(i-1)=>Q_i, Q_i=>Q_(i-1)});
    if opts.Operator === null then 
        (f-sf) // (Q_(i-1)-Q_i)
    else if opts.Operator == "Grothendieck" then
        dividedDifference((1-Q_i)*f,i)
    else error("Invalid option for Operator.")
    )


polyByDividedDifference = method(Options=>{PolyType=>"Schubert"})
polyByDividedDifference (List, PolynomialRing) := opts -> (w, Q) -> (
    n := #w;
    w0 := longestPerm(n);
    schubw0 := product(n, i->(Q_i)^(n-1-i));
    curpoly := schubw0;
    if (w0 == w) then schubw0
    else (
        v := composePerms(inverseOf(w), w0);
        redword := getOneReducedWord(v);
	if opts.PolyType == "Schubert" then
            polys := apply(reverse redword, i-> (curpoly=dividedDifference(curpoly,i); curpoly))
	else if opts.PolyType == "Grothendieck" then
	    polys = apply(reverse redword, i-> (curpoly=dividedDifference(curpoly,i, Operator=>"Grothendieck"); curpoly))
	else error ("Invalid option for PolyType.");
	polys#(#polys-1)
    )
    )

