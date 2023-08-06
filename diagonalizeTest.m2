path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");
needs ("GW-type.m2")
load ("congruenceDiagonalize.m2")
load ("congruenceDiagonalizeOverInt.m2")
load ("easyIsomorphicGW.m2")
--This file tests both diagonalize and easyIsomorphicGW and diagonalizeOverInt;

testMatrix1 = matrix(QQ,{{6/1,-1/1},{2/1,3/1}});
testMatrix2 = matrix(QQ, {{0,1/1},{1/1, 0}});
resultMatrix1 = matrix(QQ, {{1/1,0},{0,-1/1}});
resultMatrix2 = matrix(QQ, {{-1/1,0},{0,1/1}});
A=congruenceDiagonalize(testMatrix2);
E=gwClass(A);
B=gwClass(resultMatrix1);
C=gwClass(resultMatrix2);
assert(easyIsomorphicGW(E, B));
--easyIsomorphicGW works fine; returns height of solution correctly
--no solutions up to height 1, ..., i; transition matrix has height i+1;
--diagonalize works fine;

testMatrix3 = matrix(QQ, {{0,1/2},{1/2, 0}});
D=congruenceDiagonalize(testMatrix3);
F=gwClass(D);
assert(easyIsomorphicGW(F, B));

--tests congruenceDiagonalize over finite fields.
testMatrix7=matrix(GF(17), {{7, 9}, {9, 6}});
G=gwClass(congruenceDiagonalize(testMatrix7));
H=gwClass(matrix(GF(17), {{7, 0}, {0, -8}}));
assert(easyIsomorphicGW(G, H));

--------------------------------------------
--This part tests congruenceDiagonalizeOverInt;
testMatrix5=matrix(ZZ, {{0, 4},{4, 0}})
assert(congruenceDiagonalizeOverInt(testMatrix5)===matrix(ZZ, {{8, 0}, {0, -128}}));
testMatrix6=matrix(ZZ, {{6, 5},{5, 9}})
assert(congruenceDiagonalizeOverInt(testMatrix6)===matrix(ZZ, {{6, 0}, {0, 174}}));


testMatrix8=(matrix(ZZ, {{3, 6}, {6, 2}}));
assert(congruenceDiagonalizeOverInt(testMatrix8)===matrix(ZZ, {{3, 0}, {0, -90}}));
