
--Testing MatrixSchubertConstructions with Identity permutation / matrix--
TEST ///

w = {1,2,3,4};
I = matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

--isPartialASM--
assert(isPartialASM(I) == true);
assert(partialASMToASM(I) == I);

--antiDiagInit--
testIdealPerm = antiDiagInit(w);
testIdealMat = antiDiagInit(I);
assert(testIdealPerm == ideal(0_(ring testIdealPerm)));
assert(testIdealMat == ideal(0_(ring testIdealMat)));

--rankMatrix--
assert(rankTable(w) == matrix{{1,1,1,1},{1,2,2,2},{1,2,3,3},{1,2,3,4}} );
assert(rankTable(I) == matrix{{1,1,1,1},{1,2,2,2},{1,2,3,3},{1,2,3,4}} );

--rotheDiagram--
assert(rotheDiagram(w) == {} );
assert(rotheDiagram(I) == {} );

--augmentedRotheDiagram--
assert(augmentedRotheDiagram(w) == {} );
assert(augmentedRotheDiagram(I) == {} );

--essentialSet--
assert(essentialSet(w) == {} );
assert(essentialSet(I) == {} );

assert(augmentedEssentialSet(w) == {} );
assert(augmentedEssentialSet(I) == {} );

--schubDetIdeal--
testIdealPerm = schubDetIdeal(w);
testIdealMat = schubDetIdeal(I);
assert(testIdealPerm == ideal(0_(ring testIdealPerm)));
assert(testIdealMat == ideal(0_(ring testIdealMat)));

--fultonGens--
assert(fultonGens(w) == {0} );
assert(fultonGens(I) == {0} );

--diagLexInitSE--
testIdealPerm = diagLexInitSE(w);
testIdealMat = diagLexInitSE(I);
assert(testIdealPerm == monomialIdeal(0_(ring testIdealPerm)));
assert(testIdealMat == monomialIdeal(0_(ring testIdealMat)));

--diagLexInitNW--
testIdealPerm = diagLexInitNW(w);
testIdealMat = diagLexInitNW(I);
assert(testIdealPerm == monomialIdeal(0_(ring testIdealPerm)));
assert(testIdealMat == monomialIdeal(0_(ring testIdealMat)));

--diagRevLexInit--
testIdealPerm = diagRevLexInit(w);
testIdealMat = diagRevLexInit(I);
assert(testIdealPerm == monomialIdeal(0_(ring testIdealPerm)));
assert(testIdealMat == monomialIdeal(0_(ring testIdealMat)));

--subwordComplex--
--Person who wrote this needs to adjust code so that handles identity permutation as a special case.

--entrywiseMin--
--Examples in other file

--entryWiseMax--
--Examples in other file

--schubDecomposition--
testIdealPerm = schubDetIdeal(w);
testIdealMat = schubDetIdeal(I);
assert(schubDecomposition schubDetIdeal w == {{1, 2, 3, 4}})
assert(schubDecomposition schubDetIdeal I == {{1, 2, 3, 4}})

--permOverASM--
assert(permOverASM I == {{1, 2, 3, 4}})

--isIntersectionSchubIdeals--
assert(isIntersectionSchubIdeals schubDetIdeal w == true );
assert(isIntersectionSchubIdeals schubDetIdeal I == true );

--isASMIdeal--
assert(isASMIdeal schubDetIdeal w == true );
assert(isASMIdeal schubDetIdeal I == true );

--isASMUnion--
--Examples in other file

--getASM--
assert(getASM schubDetIdeal w == I );
assert(getASM schubDetIdeal I == I );

--isMinRankTable--
assert(isMinRankTable rankTable w == true );
assert(isMinRankTable rankTable I == true );

--rankTableToASM--
--Examples in other file

--rankTableFromMatrix--
--Examples in other file

--schubIntersect--
--Examples in other file

--schubAdd--
--Examples in other file

--getPermFromASM
assert(getPermFromASM getASM schubDetIdeal I == w );

--ASM--
--???

///
