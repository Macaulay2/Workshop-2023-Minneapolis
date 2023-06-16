path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");
load ("GW-type.m2")
load ("diagonalize.m2")
load ("diagonalizeOverInt.m2")
load ("easyIsomorphicGW.m2")
--This file tests both diagonalize and easyIsomorphicGW and diagonalizeOverInt;

testMatrix1 = matrix(QQ,{{6/1,-1/1},{2/1,3/1}});
testMatrix2 = matrix(QQ, {{0,1/1},{1/1, 0}});
resultMatrix1 = matrix(QQ, {{1/1,0},{0,-1/1}});
resultMatrix2 = matrix(QQ, {{-1/1,0},{0,1/1}});
A=diagonalize(testMatrix2);
E=gwClass(A);
B=gwClass(resultMatrix1);
C=gwClass(resultMatrix2);
assert(easyIsomorphicGW(E, B)===true or easyIsomorphicGW(E, C)===true);
--easyIsomorphicGW works fine; returns height of solution correctly
--no solutions up to height 1, ..., i; transition matrix has height i+1;
--diagonalize works fine;

testMatrix3 = matrix(QQ, {{0,1/2},{1/2, 0}});
D=diagonalize(testMatrix3);
F=gwClass(D);
assert(easyIsomorphicGW(F, B)===true or easyIsomorphicGW(F, C)===true);

--tests diagonalize over finite fields.
testMatrix7=matrix(ZZ/17, {{7, 9}, {9, 6}});
G=gwClass(diagonalize(testMatrix7));
H=gwClass(matrix(ZZ/17, {{7, 0}, {0, -8}}));
assert(easyIsomorphicGW(G, H)===true);

--------------------------------------------
--This part tests diagonalizeOverInt;
testMatrix5=matrix(ZZ, {{0, 4},{4, 0}})
assert(diagonalizeOverInt(testMatrix5)===matrix(ZZ, {{8, 0}, {0, -128}}));
testMatrix6=matrix(ZZ, {{6, 5},{5, 9}})
assert(diagonalizeOverInt(testMatrix6)===matrix(ZZ, {{6, 0}, {0, 174}}));


testMatrix8=(matrix(ZZ, {{3, 6}, {6, 2}}));
J=diagonalizeOverInt(testMatrix8);
print J;---This does not return the correct matrix either.
