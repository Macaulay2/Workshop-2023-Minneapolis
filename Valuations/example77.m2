needsPackage "Valuations"

R = QQ[x_1, x_2, x_3]

e_1 = x_1 + x_2 + x_3
e_2 = x_1*x_2 + x_1*x_3 + x_2*x_3
e_3 = x_1*x_2*x_3
y = (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3)

A = subring {e_1, e_2, e_3, y}

S = QQ[z_1 .. z_4] -- z_i ~ e_i, z_4 ~ y

m = map(R, S, gens A)
ker m

