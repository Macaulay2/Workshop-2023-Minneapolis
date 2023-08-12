path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");
needs "A1BrouwerDegrees.m2"
--load "globalA1Degree.m2"
--load "wittDecomp.m2"
--load "localAlgebraBasis.m2"
--load "localA1Degree.m2"
--load "isIsomorphic2.m2"
--load "diagonalize.m2"
--load "hilbertSymbol.m2"
--load "simplifyForm.m2"

T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
f1GD = globalA1Degree(f1);
f1GDmat = f1GD.matrix;
assert((wittDecomp(f1GDmat))==(3,matrix(QQ,{{}})));
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
f1LDq= localA1Degree(f1,q);
f1LDr= localA1Degree(f1,r);
f1LDsum = gwAdd(f1LDq, f1LDr);
assert(isIsomorphic2(f1LDsum, f1GD));



T2 = QQ[w];
f2 = {w^4 + w^3 - w^2 - w};
f2GD= globalA1Degree(f2);
f2GDmat = f2GD.matrix;
assert(wittDecomp(f2GDmat)==(2,matrix(QQ,{{}})));

p=ideal {w+1};
f2LDp = localA1Degree(f2, p);
f2LDpmat = f2LDp.matrix;
assert(wittDecomp(f2LDpmat)==(1,matrix(QQ,{{}})));
s=ideal{w-1};
f2LDs = localA1Degree(f2, s);
t=ideal{w};
f2LDt = localA1Degree(f2, t);
f2LDsum = gwAdd(gwAdd(f2LDp, f2LDs),f2LDt);
assert(isIsomorphic2(f2LDsum, f2GD));




ff = GF(17);
T3 = ff[y_1..y_3];
f3 = {y_1^2, y_2^2, y_3^2};
f4 = {y_2^2, y_3^2, y_1^2};
u = ideal {y_1, y_2, y_3};
f3GD = globalA1Degree(f3);
f4GD = globalA1Degree(f4);
assert(isIsomorphic2(f3GD, f4GD));
f3LDu = localA1Degree(f3, u);
f4LDu = localA1Degree(f4, u);
assert(isIsomorphic2(f3LDu, f4LDu));
assert(isIsomorphic2(f3LDu, f3GD));



T4 = CC[v];
f5 = {v^4 + v^3 - v^2 - v};
f5GD = globalA1Degree(f5);
a = ideal {v+1};
b = ideal {v};
c = ideal {v-1};
f5LDa = localA1Degree(f5,a);
f5LDb = localA1Degree(f5,b);
f5LDc = localA1Degree(f5,c);
f5LDsum = gwAdd(gwAdd(f5LDa, f5LDb),f5LDc);
assert(isIsomorphic2(f5LDsum, f5GD));
