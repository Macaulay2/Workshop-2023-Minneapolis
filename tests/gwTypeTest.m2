path = append(path, "/home/macaulay/A1-Brouwer/");
load "GW-type.m2"
print(currentDirectory());
M = matrix(QQ,{{1,0},{0,1}});
N = matrix(QQ, {{1, 2}, {3, 4}})
beta = gwClass(M);
gamma = gwClass(N);
assert(baseField(beta) === QQ)
assert(beta.matrix === M)
--Operations within GW-classes
A = gwAdd(beta, gamma);
B = gwMultiply(beta, gamma);
assert(A.matrix === matrix(QQ, {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 2}, {0, 0, 3, 4}}));
assert(B.matrix === matrix(QQ, {{1, 2, 0, 0}, {3, 4, 0, 0}, {0, 0, 1, 2}, {0, 0, 3, 4}}));
---non well-defined GW-classes
M'=matrix(ZZ, {{1, 0}, {0, 1}});
N'=matrix(QQ, {{1, 1}, {1, 1}});
theta=gwClass(M');
sigma=gwClass(N');
assert(isWellDefined(theta) === false);
assert(isWellDefined(sigma) === false);
