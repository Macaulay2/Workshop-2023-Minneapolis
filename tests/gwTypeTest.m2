load ("/home/macaulay/A1-Brouwer/GW-type.m2")
print(currentDirectory());
M = matrix(QQ,{{1,0},{0,1}});
beta = gwClass(M);
assert(baseField(beta) === QQ)
assert(beta.matrix === M)

print("got here")
