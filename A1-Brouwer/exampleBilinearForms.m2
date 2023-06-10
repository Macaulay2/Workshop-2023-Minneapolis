
b = matrix(QQ,{{1,2},{2,0}})
beta = gwClass(b)

g = matrix(QQ,{{5,6,7},{6,0,1},{7,1,0}})
gamma = gwClass(g)


-- For checking upper triangularity

b1 = matrix(QQ,{{5,6,1},{6,2,0},{1,0,0}})
b2 = matrix(GF(17), {{4,3,1,7},{3,2,5,0},{1,5,0,0},{7,0,0,0}})
g = matrix(QQ,{{5,6,7},{6,0,1},{7,1,0}})


A=matrix{{1/1, 0},{0 ,1}};
B=matrix{{0,1/1},{1 ,0 }};
C=matrix{{1,1/1},{1 ,1 }};
D=matrix{{-1/1,-1,1,1},{-1,1,1,0},{1,1,0,0},{1,0,0,0}};

print diagonalize(A)
print diagonalize(B)
print diagonalize(C)
print diagonalize(D)



-----

b1 = matrix(QQ,{{5,6,1},{6,2,0},{1,0,0}});
g = matrix(QQ,{{5,6,7},{6,0,1},{7,1,0}});
M = matrix(RR,{{0,0,1},{0,-2,0},{1,0,0}})
cong1 = matrix(QQ,{{0,0,1},{0,1,0},{1,0,0}});
cong2 = matrix(QQ,{{1,0,0},{0,-1,0},{0,0,1}});

c1 = gwClass(cong1);
c2 = gwClass(cong2);



beta = gwClass(b1);
gamma2 = gwClass(g);

easyIsomorphicGW(beta,gamma2,HeightBound=>3)


easyIsomorphicGW(c1,c2,HeightBound=>2)



M = matrix(QQ,{{0,0,1},{0,-2,0},{1,0,0}})
wittDecomp(M)



M = matrix(CC,{{0,0,1},{0,-2,0},{1,0,0}})
N = matrix(RR,{{0,1},{1,0}})

Mclass = gwClass(M)
Nclass = gwClass(N)
diagonalForm(Mclass)
diagonalForm(Nclass)



h1 = matrix(RR,{{1,0},{0,-1}});
hyp1 = gwClass(h1);
diagonalForm(hyp1)

h1an = matrix(RR,{{0,0,1},{0,2,0},{1,0,0}})

hyp1an = gwClass(h1an)
diagonalForm(hyp1an)

h2 = matrix(RR,{{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0}});
hyp2 = gwClass(h2);
diagonalForm(hyp2)



Crk3 = gwClass(matrix(CC,{{2,3},{3,5}}));
Crk3

print(simplifyForm(Crk3))


HypR = gwClass(matrix(RR,{{0,1},{1,0}}));
HypR

print(simplifyForm(HypR))





M = matrix(GF(5),{{1,0,0},{0,3,0},{0,0,4}})
beta = gwClass(M)
print(simplifyForm(beta))






T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
loc1= localA1Degree(f1,q)
loc1




T=QQ[x_1..x_4];
Endo = {x_2-x_1*x_3, x_4-x_1-x_3^2, 1-x_1*x_4, -x_2-x_3*x_4};
U= globalA1Degree(Endo);
U
diagonalize(U)
