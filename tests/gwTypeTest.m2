load "./GW-type.m2"

M = matrix(QQ,{{1,0},{0,1}});
beta = gwClass(M);
assert(baseField(beta) === QQ)
assert(beta.matrix === M)
