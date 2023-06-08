--Brianchon's thm: If a hexagon circumscribes a conic section, its diagonals intersect in a point
restart
loadPackage("Brackets",FileName=>"../Brackets.m2")

--generate the condition that the diagonals of the hexagon intersect in a single point
G = gc(a..f,3)

--the diagonals
adLine = (a*d)_G 
beLine = (b*e)_G
cfLine = (c*f)_G

--the condition
cond = (adLine ^ beLine) ^ cfLine
normalForm cond



--generate the condition that the hexagon circumscribes a conic
B = bracketRing G
R = ring B

S = QQ[gens R,q_1..q_5,z_0..z_2]
X = sub(transpose genericMatrix(R,3,6),S)

M = matrix{{1,q_1,q_2},{q_1,q_3,q_4},{q_2,q_4,q_5}}
v = matrix{{z_0},{z_1},{z_2}}
C = ((transpose v)*M*v)_(0,0)

invM = matrix table(3,3,(i,j)->(-1)^(i+j)*(det submatrix(M,delete(i,toList(0..2)),delete(j,toList(0..2)))))
dualC = ((transpose v)*invM*v)_(0,0)

L = for i from 0 to 5 list (
    dualPt = {z_0=>X_((i+1)%6,1)-X_(i,1),z_1=>X_(i,0)-X_((i+1)%6,0),z_2=>X_(i,2)-X_((i+1)%6,2)};
    sub(dualC,dualPt)
)
--crashes my computer :(
F = sub(eliminate(toList(q_1..q_5),ideal L),R)
F%(ideal B)


