# Valuations
## Setup and Motivation

### What is a valuation?

A valuation is a map \(v: R \rightarrow \Gamma\) between a ring \(R\) and an ordered semigroup \(\Gamma\) such that:

* \(v(ab) = v(a) + v(b)\)
* \(v(a+b) \ge \min\{v(a), v(b)\}\)
* (For \(k\)-algebras) \(v(\lambda a) = v(a)\)

They are fundamental for the definition of:

* Tropical varieties (polyhedral complexes that encode data about the original variety)
* Khovanskii bases (a generating set for an algebra analogous to Groebner bases)

### Let's define a valuation

Load the package:

```
needsPackage "Valuations"
```

* The \(p\)-adic valuation

The *\(p\)-adic valuation* of \(x \in \QQ\) is the highest exponent of \(p\) dividing \(x\).

```
v = padicValuation(3)
```
```
v(81)
```
```
v(1/18)
```

* The leading term valuation

The *leading term valuation* of a polynomial \(f \in K[x_1 \dots x_n]\) with respect to a monomial order is the exponent of the leading term of \(f\).

```
R = QQ[x,y,z, MonomialOrder => Lex];
v = leadTermValuation R;
```
```
f = x^3*y^2*z^2 + x^3*y*z^3 + y^10 + z^8
leadTerm f
v(f)
```

The values of the valuations are ordered:
```
g = x^4*y*z^4 + x^3 + z^4
v(g)
```
```
v(f) < v(g)
```

## Ordered \(Q\)-modules
    Valuations can return results in a custom made linearly ordered free \(Q^n\)-module type.  As an example
    ```
    R = QQ[x,y]
    f = x^2 + y
    g = y^2 + y

    val = leadTermValuation(R)
    ```

    ```
    val(f)
    ```
    ```
    val(g)
    ```

    ```
    val(f) > val(g)
    ```

    One can also create ordered \(Q^n\)-modules directly by providing the desired rank and an order to use in comparing \(Q^n\)-module elements.
    ```
    M = orderedQQn(3, {Lex})
    ```

    ```
    m_0
    ```

    ```
    m_1
    ```

    ```
    m_1 > m_2
    ```
    
## Kaveh-Manon Example

```
Load example code
```

### Valuation on \(\mathbb{C}[x_1,x_2,x_3]^{A_3}\)
* Not induced from monomial order on \(\mathbb{C}[x_1,x_2,x_3]\)
* Construction involves: tropical geometry, Khovanskii bases, invariant theory, toric geometry, ...
    
### Construction:
\[\mathbb{C}[e_1,e_2,e_3,y]\rightarrow\mathbb{C}[x_1,x_2,x_3]^{A_3}\subseteq \mathbb{C}[x_1,x_2,x_3]\]
* \(e_i\): elementary symmetric polynomial
* \(y\): Vandermonde determinant

### Step 1:
Use tropical geometry to build valuation on \(\mathbb{C}[e_1,e_2,e_3,y]\)
```
needsPackage "Valuations"
needsPackage "SubalgebraBases"
R1 = QQ[x_1, x_2, x_3]
R2 = QQ[e_1, e_2, e_3, y]
a_1 = x_1 + x_2 + x_3
a_2 = x_1*x_2 + x_1*x_3 + x_2*x_3
a_3 = x_1*x_2*x_3
w = (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3)

A = subring {a_1, a_2, a_3, w}

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
        H := gfanInitialForms(J, -1*(inCns#i), "ideal" =>true);
        K := H_1;
        if binomialIsPrime(ideal(K)) then cns#i);
    return delete(null,L))

C = primeConesOfIdeal I


coneToMatrix = coneRays -> (
    v1 := coneRays_0 + coneRays_1;
    v2 := coneRays_0 + 2*coneRays_1;
    transpose matrix {v1, v2}
    )

positivity = (f, matL) -> (
    l := transpose linealitySpace(f);
    finalScaledMats := {};
    matList := for i from 0 to #matL-1 list entries matL_i;
    
    for i from 0 to #matList-1 do (
	scaledRows = {};
	for j from 0 to #(matList_i)-1 do (
		coeff := -1*floor(min apply(#(matList_i)_j, k -> (((matList_i)_j)_k)/(flatten entries l)_k));
		scaledRows = append(scaledRows, (1/gcd(flatten entries (coeff*l + matrix{(matList_i)_j})))*(coeff*l + matrix{(matList_i)_j}));
		);
	mat := scaledRows_0;
	for i from 1 to #scaledRows-1 do mat = mat || scaledRows_i;
	finalScaledMats = append(finalScaledMats, mat);
    );
finalScaledMats
)

F = tropicalVariety(I)
M = coneToMatrix(C#1)
P = positivity(F, {M})

coneToValuation = coneRays -> (
    M := coneToMatrix(coneRays);
    scaledM := (positivity(F, {M}))/(i -> sub(i, ZZ));
    T := QQ[z_1 .. z_4, MonomialOrder=>{Weights=>((entries scaledM_0)_0), Weights=>((entries scaledM_0)_1)}];
    val := leadTermValuation(T);
    orderedM := orderedQQn(2, {Lex});
    func := (f -> (
	    valf := val(sub(f, T));
	    if valf == infinity then infinity else (
		(gens orderedM)*(scaledM_0)*(valf)
		)
	    )
	);
    valuation(func, S, orderedM)
    )

val = coneToValuation(C#0)
-- Show this off, maybe.
```

### Step 2:
Induce valuation on \(\mathbb{C}[x_1,x_2,x_3]^{A_3}\)
```
valM = (f, T, valMTwiddle) -> (
    valMfunc = (g) -> (
        R := QQ[x_1, x_2, x_3, e_1, e_2, e_3, y, MonomialOrder => Eliminate 3];
        I := ideal{x_1 + x_2 + x_3 - e_1, x_1*x_2 + x_1*x_3 + x_2*x_3 - e_2, x_1*x_2*x_3 - e_3, (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3) - y};
        S := valMTwiddle#"domain";
	m := map(S, R, matrix{{0,0,0}} | vars S);
        gTwiddle := m (sub(g, R) % I);
        maxTwiddle := gTwiddle % ideal(sub(f, S));
        valMTwiddle(maxTwiddle)
    );
    valuation(valMfunc, T, valMTwiddle#"codomain")
    )
```

### Example:
```
f = e_1^2*e_2^2 - 4*e_2^3 - 4*e_3*e_1^3 + 18*e_1*e_2*e_3 - 27*e_3^2 - y^2;
W = QQ[x_1, x_2, x_3]
finalValuation = valM(f, W, val)
use W
finalValuation(x_1 + x_2 + x_3)
use W
finalValuation((x_1^2 - x_2^2)*(x_1^2 - x_3^2)*(x_2^2 - x_3^2))
use W
finalValuation(0_W)
use W
finalValuation(x_1)
use W
finalValuation(x_1*x_2*x_3)
-- val = coneToValuation(C#0) can also do this for C#1 and C#2
```


## Future Directions
1. General development
2. Use valuations in other M2 pacakges:
   - [Tropical](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/Tropical/html/index.html) --  the main M2 package for tropical computations
   - [FormalGroupLaws](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/FormalGroupLaws/html/index.html) -- commutative formal group laws
   - [SubalgebraBases](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.18/share/doc/Macaulay2/SubalgebraBases/html/index.html) -- A package for finding canonical subalgebra bases (Sagbi bases)
3. Puiseux series
4. Any suggestions?
