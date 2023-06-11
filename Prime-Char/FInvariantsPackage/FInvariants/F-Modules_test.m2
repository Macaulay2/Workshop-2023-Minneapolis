-------------------
-- 24 hours of Cohomology, Example 21.31

TEST ///
R = ZZ/5[x_1..x_6];
I = minors( 2, matrix {{x_1,x_2,x_3},{x_4,x_5,x_6}} );
J = minors( 2, matrix {{x_1,x_2,x_3},{x_4,x_5,x_6}} );
assert( cohomDim I == 2 )
assert( cohomDim( J, Strategy => Filter ) == 2 )
///
  
-------------------
-- Theorem 3.2 in https://arxiv.org/abs/1503.06184v1

TEST ///
n = (2,2,3);
d = #n;
variables = fold( apply(d, i -> x_(i+1,0)..x_(i+1,n#i)), (i,j) -> i|j );
R = ZZ/2[variables];
B = apply(d, i -> transpose matrix toList apply(n#i, j-> { x_(i+1,j), x_(i+1,j+1) }) );
B = fold( B, (a,b) -> a | b );
J = minors( 2, B );
assert( cohomDim J == sum( toList n ) - 1 )
///

-------------------
-- Result mentioned in intro of E.E. Witt / Advances in Mathematics 231 (2012) 1998–2012 

TEST ///
r = 2;
s = 6;
R = ZZ/3[x_(1,1)..x_(r,s)];
X = matrix toList apply(1..r, i -> toList apply(1..s, j -> x_(i,j) ) );
I = minors( min(r,s), X );
assert( all( 0..(s-r), i -> localCohomology( i, I, R ) == 0 ) and cohomDim I == s - r + 1 )
///

-----------------
-- Example 4.2 in https://arxiv.org/pdf/1509.01519.pdf

TEST ///
R = ZZ/2[x_0,x_1,x_2,y_0,y_1,y_2,z_0,z_1,z_2];
S = R[t];
F1 = x_0+x_1*t+x_2*t^2;
F2 = y_0+y_1*t+y_2*t^2;
F3 = z_0+z_1*t+z_2*t^2;
r1 = resultant( F2, F3, t );
r2 = resultant( F1, F3, t );
r3 = resultant( F1, F2, t );
r4 = resultant( F1 + F2, F3, t );
I = ideal( r1, r2, r3, r4 );
assert( localCohomology( 4, I, R ) == 0 )
///

------------------
-- Examples from Katzman--Sharp, TAMS 375 No. 9, 2022, pp. 6621--6651

TEST ///
R = ZZ/2[x_1..x_4];
I = ideal(x_1*x_2, x_2*x_3, x_3*x_4, x_4*x_1); -- Example 6.7(ii)
assert( lyubeznikTable( I, R ) == matrix { {0,1,0}, {0,0,0}, {0,0,2} } )
///

TEST ///
R = ZZ/2[x_1..x_5];
I = intersect(ideal(x_1..x_3), ideal(x_3..x_5), ideal(x_1..x_4)); -- Example 6.7(iii)
assert( lyubeznikTable( I, R ) == matrix { {0,1,0}, {0,0,0}, {0,0,2} } )
///

TEST ///
R = ZZ/2[x_1..x_6];
I = ideal(x_1*x_2*x_3, x_1*x_2*x_4, x_1*x_3*x_5, x_2*x_4*x_5, x_3*x_4*x_5, x_2*x_3*x_6, x_1*x_4*x_6, x_3*x_4*x_6, x_1*x_5*x_6, x_2*x_5*x_6); -- Example 6.7(i)
assert( lyubeznikTable( I, R ) == matrix { {0,0,1,0}, {0,0,0,0}, {0,0,0,1}, {0,0,0,1} } )
J = intersect(ideal(x_1,x_2), ideal(x_3,x_4), ideal(x_5,x_6)); -- Exampl 6.7(iv)
assert( lyubeznikTable( J, R ) == matrix { {0,0,1,0,0}, {0,0,0,0,0}, {0,0,0,3,0}, {0,0,0,0,0}, {0,0,0,0,3} } )
K = ideal( x_1*x_2, x_3*x_4, x_5*x_6, x_1*x_3*x_5, x_2*x_4*x_6 ); -- Example 6.7(vi)
assert( lyubeznikTable( K, R ) == matrix { {0,0,1,0}, {0,0,0,0}, {0,0,0,1}, {0,0,0,1} } )
///

TEST ///
R = ZZ/2[x_1..x_7];
I = intersect(ideal(x_1,x_2), ideal(x_3,x_4), ideal(x_5..x_7)); -- Example 6.6
assert( lyubeznikTable( I, R ) == matrix { {0,0,1,0,0,0}, {0,0,0,0,0,0}, {0,0,0,2,0,0}, {0,0,0,0,1,0}, {0,0,0,0,1,0}, {0,0,0,0,0,2}} )
J = ideal( x_1*x_2, x_2*x_3, x_3*x_4, x_4*x_5, x_5*x_6, x_6*x_7, x_7*x_1 ); -- Example 6.7(v)
assert( lyubeznikTable( J, R ) == matrix { {0,0,1,0}, {0,0,0,0}, {0,0,0,1}, {0,0,0,1} } )
///

