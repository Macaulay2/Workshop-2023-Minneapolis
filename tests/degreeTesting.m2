path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");
load "globalA1Degree.m2"
load "wittDecomp.m2"
load "localAlgebraBasis.m2"
load "localA1Degree.m2"

T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
assert((wittDecomp(globalA1Degree(f1)))==(3,matrix(QQ,{{}})));
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
loc1= localA1Degree(f1,q);
loc2= localA1Degree(f1,r);
assert((wittDecomp(loc1++loc2))==(3,matrix(QQ,{{}})));

T2 = QQ[w];
f2 = {w^4 + w^3 - w^2 - w};
assert( wittDecomp(globalA1Degree(f2))==(2,matrix(QQ,{{}})));
p=ideal {w+1/1};
assert(wittDecomp(localA1Degree(f2,p))==(1,matrix(QQ,{{}})));
s=ideal{w-1};
t=ideal{w};
assert(wittDecomp(localA1Degree(f2,t)++localA1Degree(f2,s))==(1,matrix(QQ,{{}})));

ff = ZZ/17;
T3 = ff[y_1..y_3];
f3 = {y_1^2, y_2^2, y_3^2};
f4 = {y_2^2, y_3^2, y_1^2};
u = ideal {y_1, y_2, y_3};
assert(globalA1Degree(f3) == localA1Degree(f3, u));
assert(globalA1Degree(f4) == localA1Degree(f4, u));
assert(globalA1Degree(f3) == globalA1Degree(f4));
assert(globalA1Degree(f3) == localA1Degree(f4, u));
assert(globalA1Degree(f4) == localA1Degree(f3, u));