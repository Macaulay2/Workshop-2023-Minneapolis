load ("/home/macaulay/A1-Brouwer/GW-type.m2")
load ("/home/macaulay/A1-Brouwer/diagonalize.m2")
load ("/home/macaulay/A1-Brouwer/easyIsomorphicGW.m2")

testMatrix1 = matrix(QQ,{{6/1,-1/1},{2/1,3/1}});
testMatrix2 = matrix(QQ, {{0,1/1},{1/1, 0}});
resultMatrix1 = matrix(QQ, {{1/1,0},{0,-1/1}});
resultMatrix2 = matrix(QQ, {{-1/1,0},{0,1/1}});
A=diagonalize(testMatrix2);
E=gwClass(A);
B=gwClass(resultMatrix1);
C=gwClass(resultMatrix2);
