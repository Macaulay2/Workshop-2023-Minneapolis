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
time cohomDim J == 2
time cohomDim( J, Strategy => Filter ) == 2

-------------------
-- Theorem 3.2 in https://arxiv.org/abs/1503.06184v1

p = 2;
n = 2;
d = 3;
R = ZZ/p[x_(1,0)..x_(d,n)];
B = apply(1..d, i -> transpose matrix toList apply(n, j-> { x_(i,j), x_(i,j+1) }) );
B = fold( B, (a,b) -> a | b )
J = minors( 2, B );

cohomDim J == n*d - 1

time cohomDim( J, Strategy => Filter ) == n*d - 1 -- Takes too long!

p = 2;
n = (2,2,3);
d = #n;
N = max(n)
variables = fold( apply(d, i -> x_(i+1,0)..x_(i+1,n#i)), (i,j) -> i|j )
R = ZZ/p[variables]

B = apply(d, i -> transpose matrix toList apply(n#i, j-> { x_(i+1,j), x_(i+1,j+1) }) )
B = fold( B, (a,b) -> a | b )

J = minors( 2, B );

cohomDim J == sum( toList n ) - 1

-------------------
-- Result mentioned intro of in E.E. Witt / Advances in Mathematics 231 (2012) 1998–2012 

p = 2;
r = 3;
s = 6;
R = ZZ/p[x_(1,1)..x_(r,s)]
X = matrix toList apply(1..r, i -> toList apply(1..s, j -> x_(i,j) ) );
I = minors( min(r,s), X );
time cohomDim I == s - r + 1 and all( 0..(s-r), i -> localCohomology( i, I, R ) == 0 )
-- 11s
time cohomDim( I, Strategy => Filter ) == s - r + 1 and all( 0..(s-r), i -> localCohomology( i, I, R, Strategy => Filter ) == 0 )  -- takes too long

-----------------
-- Example 4.2 in https://arxiv.org/pdf/1509.01519.pdf

R = ZZ/2[x_0,x_1,x_2,y_0,y_1,y_2,z_0,z_1,z_2];
S = R[t];
F1 = x_0+x_1*t+x_2*t^2;
F2 = y_0+y_1*t+y_2*t^2;
F3 = z_0+z_1*t+z_2*t^2;
r1 = resultant( F2, F3, t );
r2 = resultant( F1, F3, t );
r3 = resultant( F1, F2, t );
r4 = resultant( F1 + F2, F3, t );
I = ideal( r1, r2, r3, r4 )

localCohomology( 4, I, R ) == 0 

f = inducedMap(R^1/I, R^1/(frobenius I) )

E = Ext^4( f, R^1 )


