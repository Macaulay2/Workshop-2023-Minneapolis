restart
needsPackage "Valuations"

R = QQ[x_1, x_2, x_3]

e_1 = x_1 + x_2 + x_3
e_2 = x_1*x_2 + x_1*x_3 + x_2*x_3
e_3 = x_1*x_2*x_3
y = (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3)

A = subring {e_1, e_2, e_3, y}

S = QQ[z_1 .. z_4] -- z_i ~ e_i, z_4 ~ y

m = map(R, S, gens A)
I = ker m

needsPackage "Tropical"
needsPackage "gfanInterface"
needsPackage "Binomials"

primeConesOfIdeal = I -> (F:=tropicalVariety(I, IsHomogeneous=>true,Prime=>true);
    r:=rays(F);
    c:=maxCones(F);
    cns := for i in c list(r_i);
    inCns := for c in cns list (flatten entries( c * transpose matrix{toList(numColumns(c) : 1)}));
    L:= for i from 0 to #cns-1 list (J = gfanBuchberger(I, "w" => -1*(inCns#i));
        H = gfanInitialForms(J, -1*(inCns#i), "ideal" =>true);
        K = H_1;
        if binomialIsPrime(ideal(K)) then K);
    return delete(null,L))

primeConesOfIdeal I