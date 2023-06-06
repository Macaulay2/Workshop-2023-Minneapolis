--Brianchon's thm: If a hexagon circumscribes a conic section, its diagonals intersect in a point
restart
loadPackage("Brackets",FileName=>"../Brackets.m2")
G = gc(a..f,3)
B = bracketRing G
R = ring B

S = QQ[(gens R),q_1..q_5,t,z_0..z_2]
X = sub(transpose genericMatrix(R,3,6),S)

coeffs = transpose (matrix{{1}|toList(q_1..q_5)})
C = (basis(2,S,Variables=>toList(z_0..z_2))*coeffs)_(0,0)

DList = for i from 0 to 5 list (
    linPara = apply(3,j->z_j=>t*X_(i,j)+(1-t)*X_((i+1)%6,j));
    Csub = sub(C,linPara);
    discriminant(Csub,t)
)
myIdeal = eliminate(toList(q_1..q_5),ideal DList)


