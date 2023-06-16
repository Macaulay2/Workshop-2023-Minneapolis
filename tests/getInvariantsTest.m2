load "./diagonalize.m2"
loadPackage "RationalPoints2"
load "GW-type.m2"

-- test code and assertions
needsPackage "getInvariants"

alpha = gwClass(matrix(QQ,{{1,2,2},{2,2,0},{2,0,0}}))
beta = gwClass(matrix(QQ,{{2,4,4},{4,4,0},{4,0,0}}))

getInvariants(alpha)
getInvariants(beta)

assert(getInvariants(alpha) == {3, {2, 0, 1}, -8})
assert(isIdenticalDiscriminant(alpha,beta) == false)
