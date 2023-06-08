-- Pascal's Theorem

restart
debug needsPackage "Brackets"

-*
Under what conditions do 6 points in the projective plane lie on a single quadric?

Choose six projective points a, b, c, d, e, f
*-

G = gc(a .. f, 3)
gens G
B = bracketRing G;
X = matrix B;
C = fold(apply(0..5, i-> basis(2, ring X, Variables => (entries X)#i)), (a,b) -> a||b); -- Matrix with rows corresponding to six point on a quadric
D = det C; -- D = 0 if and only if the six points lie on a single conic
q1 = toBracketPolynomial(D, B)

abLine = (a * b)_G -- Line joining a and b
afLine = (a * f)_G -- Line joining a and f
edLine = (e * d)_G -- Line joining e and d
efLine = (e * f)_G -- Line joining e and f
cdLine = (c * d)_G -- Line joining c and d
bcLine = (b * c)_G -- Line joining b and c

p1 = abLine ^ edLine -- Intersection point of lines joining a, b and e, d
p2 = afLine ^ cdLine -- Intersection point of lines joining a, f and c, d
p3 = bcLine ^ efLine -- Intersection point of lines joining b, c and e, f

q2 = p1 * p2 * p3 -- Span of p1, p2, p3. q = 0 if the points are collinear.

normalForm q2 === (-1) * q1 -- True! So, a,b,c,d,e,f lie on a single quadric if and only if p1, p2, p3 are collinear.

-*
To document:

  gc
  gens
  bracketRing  DONE
  matrix 
  toBracketPolynomial
  * 
  _
  ^
  normalForm
*-

