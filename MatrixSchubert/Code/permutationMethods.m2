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
