R=ZZ/101[x_1..x_12]
I=ideal(gens R)
cohomDim I

-------------------

p = 2
R = ZZ/p[x_1..x_7];
I = intersect(ideal(x_1, x_2), ideal(x_3, x_4), ideal(x_5, x_6, x_7));
#(trim I)_*
cohomDim I

J = minors( 2, matrix {{0,x_1,x_2,x_3},{x_4,x_5,x_6,x_7}} )
cohomDim J

-------------------
-- 24 hours of Cohomology, Example 21.31

p = 5
R = ZZ/p[x_1..x_6]
J = minors( 2, matrix {{x_1,x_2,x_3},{x_4,x_5,x_6}} )
cohomDim J == 2

-------------------
-- Theorem 3.2 in https://arxiv.org/abs/1503.06184v1

p = 2;
n = 2;
d = 3;
R = ZZ/p[x_(1,0)..x_(d,n)];
B = apply(1..d, i -> transpose matrix toList apply(n, j-> { x_(i,j), x_(i,j+1) }) );
B = fold( B, (a,b) -> a | b );
J = minors( 2, B );
cohomDim J == n*d - 1
