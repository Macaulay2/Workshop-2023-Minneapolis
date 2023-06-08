restart
needsPackage "Valuations"
needsPackage "SubalgebraBases"


R = QQ[x_1, x_2, x_3, e_1, e_2, e_3, y, MonomialOrder => Eliminate 3]
I = ideal{x_1 + x_2 + x_3 - e_1, x_1*x_2 + x_1*x_3 + x_2*x_3 - e_2, x_1*x_2*x_3 - e_3, (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3) - y}
f = e_1^2*e_2^2 - 4*e_2^3 - 4*e_3*e_1^3 + 18*e_1*e_2*e_3 - 27*e_3^2 - y^2

valM = (f, g, valMTwiddle) -> (
    R : = QQ[x_1, x_2, x_3, e_1, e_2, e_3, y, MonomialOrder => Eliminate 3];
    I := ideal{x_1 + x_2 + x_3 - e_1, x_1*x_2 + x_1*x_3 + x_2*x_3 - e_2, x_1*x_2*x_3 - e_3, (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3) - y};
    S := valMTwiddle#"domain"
    gTwiddle = sub(sub(g, R) % I, S);
    maxTwiddle = gTwiddle % ideal(sub(f, S));
    valMTwiddle(maxTwiddle)
    )

-----


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

primeConesOfIdeal = I -> (
    F:=tropicalVariety(I, IsHomogeneous=>true,Prime=>true);
    r:=rays(F);
    c:=maxCones(F);
    cns := for i in c list(r_i);
    inCns := for c in cns list (flatten entries( c * transpose matrix{toList(numColumns(c) : 1)}));
    L:= for i from 0 to #cns-1 list (J = gfanBuchberger(I, "w" => -1*(inCns#i));
        H = gfanInitialForms(J, -1*(inCns#i), "ideal" =>true);
        K = H_1;
        if binomialIsPrime(ideal(K)) then cns#i);
    return delete(null,L))

primeConesOfIdeal I

coneToMatrix = coneRays -> (
    v1 = coneRays_0 + coneRays_1;
    v2 = coneRays_0 + 2*coneRays_1;
    transpose matrix {v1, v2}
    )
