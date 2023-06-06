-- pascal's theorem
restart
loadPackage("Brackets",FileName=>"../Brackets.m2")
G = gc(a .. f, 3)

abLine = (a * b)_G
afLine = (a * f)_G
edLine = (e * d)_G
efLine = (e * f)_G
cdLine = (c * d)_G
bcLine = (b * c)_G

p1 = abLine ^ edLine
p2 = afLine ^ cdLine
p3 = bcLine ^ efLine

q = p1 * p2 * p3
-- a*b*c is a top-exterior power so it's just the determinant [abc]
-- Improvement: interpret a*b*c as [abc] here.

-- peek under the hood
t = flatten ((x -> first entries (coefficients x)_1) \ terms q#RingElement)
R = ring G#bracketRing
use R
tR = t / (x -> sub(x, R))
V = vars R

q' = tR_0 * V_(0,37) + tR_1 * V_(0,34) + tR_2 * V_(0,31) + tR_3 * V_(0,25)
I = ideal G#bracketRing

normForm = q' % I
normForm_(bracketRing ring q) -- this is what we want automatically




--thomas trying math
B = bracketRing G
X = matrix B

C = fold(apply(0..5,i->basis(2,ring X,Variables => (entries X)#i)),(a,b)->a||b)
D = det C
toBracketPolynomial(D)
I = B#ideal
D % I --looks pretty
