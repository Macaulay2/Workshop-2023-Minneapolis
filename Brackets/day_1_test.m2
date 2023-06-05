restart

needsPackage "Brackets"

B = bracketRing(6, 3)
T = [1 4 5]_B * [1 5 6]_B * [2 3 4]_B
n = normalForm T 


G = gc(a..f,3) -- Grassmann-Cayley algebra for 6 points in P^2
abLine = (a * b)_G -- line spanned by a and b
deLine = (d * e)_G -- line spanned by d and e
bcLine = (b * c)_G -- line spanned by b and c
efLine = (e * f)_G -- line spanned by e and f
acLine = (a * c)_G -- line spanned by a and c
dfLine = (d * f)_G -- line spanned by d and f
pt1 = abLine ^ deLine -- intersection of ab and de
pt2 = bcLine ^ efLine -- intersection of bc and ef
pt3 = acLine ^ dfLine -- intersection of ac and df
linePerspective = pt1 * pt2 * pt3 -- Condition that the pts p1, p2, p3 are collinear
adLine = (a * d)_G -- line spanned by a and d
beLine = (b * e)_G -- line spanned by b and e
cfLine = (c * f)_G -- line spanned by c and f
pointPerspective =  adLine ^ beLine ^ cfLine -- Condition that the 3 lines above meet.
-*
"pointPerspective" and "linePerspective" are two degree-0 elements of the Grassmann Cayley
algebra, which We identify with elements of the bracket ring B_(2,6).
The representatives of these elements do not share any common factors.
But, applying the straightening algorithm below produces two normal forms such that
 [abc] * [def] * nf(linePerspective) = 2 * nf(pointPerspective)
So, if a,b,c are not collinear and d,e,f are not collinear,
line and point perspective are the same.
*-
factor linePerspective
factor pointPerspective
(n1, n2) = (normalForm pointPerspective, normalForm linePerspective);
netList factor n1
netList factor n2
