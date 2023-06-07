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

TEST ///
-- isMinRankTable
T1 = matrix {{0,1,1},{1,1,2},{1,2,3}}
T2 = matrix {{1,1,1,1,1},{1,2,2,2,3},{1,2,2,2,3},{1,2,2,3,3},{1,2,3,3,3}}
F1 = matrix {{1,0,1,0},{0,1,0,-1},{2,2,0,0},{3,5,8,0}}
F2 = matrix {{1,1,1,1,1},{1,2,2,2,0},{1,2,2,2,3},{1,2,2,3,3},{1,2,3,3,3}}

assert(isMinRankTable(T1))
assert(isMinRankTable(T2))
assert(not isMinRankTable(F1))
assert(not isMinRankTable(F2))
///
