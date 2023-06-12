-------------------
-- 24 hours of Cohomology, Example 21.31

TEST ///
R = ZZ/5[x_1..x_6];
I = minors( 2, matrix {{x_1,x_2,x_3},{x_4,x_5,x_6}} );
J = minors( 2, matrix {{x_1,x_2,x_3},{x_4,x_5,x_6}} ); 
-- need a "different" ideal, or else the second call of cohomDim will just use cached result
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
///

TEST ///
R = ZZ/2[x_1..x_6];
J = intersect(ideal(x_1,x_2), ideal(x_3,x_4), ideal(x_5,x_6)); -- Exampl 6.7(iv)
assert( lyubeznikTable( J, R ) == matrix { {0,0,1,0,0}, {0,0,0,0,0}, {0,0,0,3,0}, {0,0,0,0,0}, {0,0,0,0,3} } )
///

TEST ///
R = ZZ/2[x_1..x_6];
K = ideal( x_1*x_2, x_3*x_4, x_5*x_6, x_1*x_3*x_5, x_2*x_4*x_6 ); -- Example 6.7(vi)
assert( lyubeznikTable( K, R ) == matrix { {0,0,1,0}, {0,0,0,0}, {0,0,0,1}, {0,0,0,1} } )
///

TEST ///
R = ZZ/2[x_1..x_7];
I = intersect(ideal(x_1,x_2), ideal(x_3,x_4), ideal(x_5..x_7)); -- Example 6.6
assert( lyubeznikTable( I, R ) == matrix { {0,0,1,0,0,0}, {0,0,0,0,0,0}, {0,0,0,2,0,0}, {0,0,0,0,1,0}, {0,0,0,0,1,0}, {0,0,0,0,0,2}} )
///

TEST ///
R = ZZ/2[x_1..x_7];
J = ideal( x_1*x_2, x_2*x_3, x_3*x_4, x_4*x_5, x_5*x_6, x_6*x_7, x_7*x_1 ); -- Example 6.7(v)
assert( lyubeznikTable( J, R ) == matrix { {0,0,1,0}, {0,0,0,0}, {0,0,0,1}, {0,0,0,1} } )
///

