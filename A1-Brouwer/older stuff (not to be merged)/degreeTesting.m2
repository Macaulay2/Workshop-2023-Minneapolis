load "globalA1Degree.m2"
load "localA1Degree.m2"
load "WittDecomp.m2"

T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
assert((WittDecomp(globalA1Degree(f1)))==(3,matrix(QQ,{{}})));
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
loc1= localA1Degree(f1,q);
loc2= localA1Degree(f1,r);
assert((WittDecomp(loc1++loc2))==(3,matrix(QQ,{{}})));

T2 = QQ[w];
f2 = {w^4 + w^3 - w^2 - w};
assert( WittDecomp(globalA1Degree(f2))==(2,matrix(QQ,{{}})));
p=ideal {w+1/1};
assert(WittDecomp(localA1Degree(f2,p))==(1,matrix(QQ,{{}})));
s=ideal{w-1};
t=ideal{w};
assert(WittDecomp(localA1Degree(f2,t)++localA1Degree(f2,s))==(1,matrix(QQ,{{}})));

