load ("/home/macaulay/A1-Brouwer/GW-type.m2")
load ("/home/macaulay/A1-Brouwer/diagonalize.m2")
load ("/home/macaulay/A1-Brouwer/easyIsomorphicGW.m2")
--This file tests both diagonalize and easyIsomorphicGW;

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

testMatrix4=matrix(RR, {{0,1/2},{1/2, 0}});