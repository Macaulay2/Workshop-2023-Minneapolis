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
assert(sort essentialSet(A) == {(1,3),(2,1),(3,2)})
///

--test for schubDetIdeal
--TODO: make more complicated tests
TEST ///
--Example 15.4 from Miller-Sturmfels
I = schubDetIdeal({1,2,3});
assert(I == ideal(0_(ring I)));

I = schubDetIdeal({2,1,3});
assert(I == ideal((ring I)_0));

I = schubDetIdeal({2,3,1});
assert(I == ideal((ring I)_0, (ring I)_3));

I = schubDetIdeal({3,2,1});
assert(I == ideal((ring I)_0, (ring I)_1, (ring I)_3));

I = schubDetIdeal({1,3,2});
assert(I == ideal((ring I)_0 * (ring I)_4 - (ring I)_1 * (ring I)_3));
///

-*
--NOTE: it would be nice to keep this test, but it would also be nice to not 
--      have to export it as well
TEST ///
--indexOfVariable
R = QQ[x_1..x_5]
assert(indexOfVariable x_1 == 1)

R = QQ[x_1..x_3,y_3..y_5]
assert(indexOfVariable y_4 == 4)

R = QQ[x_(1,1)..x_(4,4)]
assert(indexOfVariable x_(2,3) == (2,3))

R = QQ[x_{1,1}..x_{4,4}]
assert(indexOfVariable x_{2,3} == {2,3})
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
--isCDG
assert(isCDG({5,4,3,2,1}));

--isCartwrightSturmfels
assert(isCartwrightSturmfels({5,4,3,2,1}));
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
-- rotheDiagram

assert(partialASMToASM matrix{{0,0,1,0},{1,0,-1,0},{0,0,0,0}} == matrix{{0,0,1,0,0,0},{1,0,-1,0,1,0},{0,0,0,0,0,1},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0}})
///

TEST ///
-- rotheDiagram

assert(sort rotheDiagram matrix{{0,0,1,0,0},{1,0,-1,0,0},{0,1,0,0,0},{0,0,1,0,0}} == sort {(1,1),(1,2),(2,3),(2,4),(2,5)})
assert(sort rotheDiagram {2,6,5,1,4,3} == sort {(1,1),(2,1),(2,3),(2,4),(2,5),(3,1),(3,3),(3,4),(5,3)})
///

TEST ///
-- augmentedRotheDiagram 

assert(sort augmentedRotheDiagram {2,1,5,4,3} == sort {((1,1),0), ((3,3),2),((3,4),2), ((4,3),2)})
assert(sort augmentedRotheDiagram matrix{{0,1,0},{1,-1,1},{0,1,0}} == sort{((1,1),0), ((2,2),1)})
assert (sort augmentedRotheDiagram matrix {{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}} == sort {((1,1),0),((1,2),0),((4,3),2),((3,3),2)})
///


TEST ///
-- essentialBoxes 

assert(essentialSet({2,1,6,3,5,4 })== {(1, 1), (3, 5), (5, 4)})
assert(essentialSet matrix {{0,1,0,0,0,0},{1,0,0,0,0,0},{0,0,0,0,0,1},{0,0,1,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0}} == {(1, 1), (3, 5), (5, 4)})

assert(essentialSet({1}) == {})
assert(essentialSet(matrix {{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}}) == {(1,2),(4,3)}) -- previously broken example
///


TEST ///
-- schubCodim
L = {
    {1},
    {2,1},
    matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}},
    matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}}
}
assert all (L, i -> schubCodim i == codim schubDetIdeal i)
///

TEST ///
--schubReg
L = {
    {2,9,3,4,1,7,5,6,8}, -- example 1.2 in PSW
    {2,1},
    matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}},
    matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}},
    matrix{{0,1,0,0},{1,0,0,0},{0,0,0,1},{0,0,1,0}}
}
expected = {5,0,1,4,2}

assert all (#L, i -> schubReg L#i == expected#i)
///

TEST ///
-- isMinRankTable
T1 = matrix {{0,1,1},{1,1,2},{1,2,3}}
T2 = matrix {{1,1,1,1,1},{1,2,2,2,2},{1,2,2,2,3},{1,2,2,3,3},{1,2,3,3,3}}
F1 = matrix {{1,0,1,0},{0,1,0,-1},{2,2,0,0},{3,5,8,0}}
F2 = matrix {{1,1,1,1,1},{1,2,2,2,0},{1,2,2,2,3},{1,2,2,3,3},{1,2,3,3,3}}

assert(isMinRankTable(T1))
assert(isMinRankTable(T2))
assert(not isMinRankTable(F1))
assert(not isMinRankTable(F2))
///

TEST///
-- rankTableToASM
Ar = matrix {{0,0,1,1},{0,1,1,2},{1,2,2,3},{1,2,3,4}}
A = matrix {{0,0,1,0},{0,1,-1,1},{1,0,0,0},{0,0,1,0}}
assert(rankTableToASM(Ar) == A)

Br = matrix {{0,0,1,1,1},{1,1,1,2,2},{1,2,2,3,3},{1,2,3,4,4},{1,2,3,4,5}}
B = matrix {{0,0,1,0,0},{1,0,-1,1,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,0,1}}
assert(rankTableToASM(Br) == B)
///

TEST///
-- rankTableFromMatrix
Am = matrix {{1,0,0},{0,23,24},{23,24,25}}
A = matrix {{0,0,0},{0,1,1},{1,2,2}}
assert(rankTableFromMatrix(Am) == A)

///

TEST ///
-- getPermFromASM
L = {
    matrix{{1,0,0,0}},
    matrix{{1,0},{0,0}},
    matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}},
    matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}},
    matrix{{0,1,0,0},{1,0,0,0},{0,0,1,0},{0,0,0,1}}
}
Lc = {{},{},{},{},{2,1,3,4}}
assert all (#L, i -> getPermFromASM L_i == Lc_i)
///
