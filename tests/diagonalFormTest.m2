path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");
load "diagonalForm.m2"
--load "GW-type.m2"

M1=matrix(RR, {{0, 1}, {1, 0}});
G1=gwClass(M1);
M2=diagonalForm(G1);
assert(M2.matrix===matrix(RR, {{1, 0}, {0, -1}}));
------------------------
M3=matrix(CC, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
assert(M4.matrix===matrix(CC, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
------------------------

M5=matrix(RR, {{1, 2, 3, 4}, {2, 4, 5, 16}, {3, 5, 7, 8}, {4, 16, 8, 19}});
G3=gwClass(M5);
M6=diagonalForm(G3);
A=matrix(RR, {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}});
assert(M6.matrix===A);
-------testings for R and C were successful;

M3=matrix(QQ, {{1, 2, 3}, {2, 4, 5}, {3, 5, 7}});
G2=gwClass(M3);
M4=diagonalForm(G2);
print M4.matrix;
---------Testing for Q was not successful.
