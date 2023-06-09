# Valuations
## Introduction and Motivation

### What is a valuation?

Valuation: map \(v: R \rightarrow \Gamma\)

\(R\): ring

\(\Gamma\): ordered semigroup

* \(v(ab) = v(a) + v(b)\)
* \(v(a+b) \ge \min\{v(a), v(b)\}\)
* \(v(x)=\infty\) iff \(x=0\)
* (For \(k\)-algebras) \(v(\lambda a) = v(a)\)

Fundamental for:

* Tropical varieties (polyhedral shadows of varieties)
* Khovanskii bases (algebra analogues of Groebner bases)

### Let's define a valuation

```
restart
needsPackage "Valuations"
viewHelp "Valuations"
```

*\(p\)-adic valuation* of \(x \in \QQ\): power of \(p\) in \(x\)

```
v = padicValuation(3)
```
```
v(81)
```
```
v(1/18)
```

*leading term valuation* of a polynomial: exponent of leading term


```
R = QQ[x,y,z, MonomialOrder => Lex];
v = leadTermValuation R;
netList pairs v
```
```
f = x^3*y^2*z^2 + x^3*y*z^3 + y^10 + z^8
leadTerm f
v(f)
```

## Ordered \(\mathbb{Q}\)-modules

Valuations can return in custom linearly ordered \(\mathbb{Q}^n\)-module type

```
R = QQ[x,y];
f = x^2 + y;
g = y^2 + y;
val = leadTermValuation(R);
```

```
val(f)
val(g)
```

```
val(f) > val(g)
```

Creation from rank and monomial ordering
```
M = orderedQQn(3, {Lex})
```

```
M_0
M_1
```

```
M_0 < M_1
```
    
## Kaveh-Manon Example

### Valuation on \(\mathbb{C}[x_1,x_2,x_3]^{A_3}\)
* Not induced from monomial order on \(\mathbb{C}[x_1,x_2,x_3]\)
* Construction involves: tropical geometry, Khovanskii bases, invariant theory, ...
    
### Construction:
\[\mathbb{C}[e_1,e_2,e_3,y]\rightarrow\mathbb{C}[x_1,x_2,x_3]^{A_3}\subseteq \mathbb{C}[x_1,x_2,x_3]\]
* \(e_i\): elementary symmetric polynomial
* \(y\): Vandermonde determinant

```
load "example77.m2"
S = QQ[e_1, e_2, e_3, y];
R = QQ[x_1, x_2, x_3];
-- SubalgebraBases package
A = subring {
    x_1 + x_2 + x_3,
    x_1*x_2 + x_1*x_3 + x_2*x_3,
    x_1*x_2*x_3,
    (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3)}; 
presMap = map(R, S, gens A);
I = ker presMap
```

### Step 1:
Use tropical geometry to build valuation on \(\mathbb{C}[e_1,e_2,e_3,y]\)
```
--- Tropical and BinomialIdeals packages
C = primeConesOfIdeal I
flatten (C/coneToMatrix/(i -> 
     positivity(tropicalVariety I, {i})))
```
```
v0 = coneToValuation(C#0, I);
v1 = coneToValuation(C#1, I);
v2 = coneToValuation(C#2, I);
```
```
use S;
f = e_1^2 + e_2*e_3 - y^3
v0(f) -- lead term from e_1^2
v1(f) -- lead term from y^3
v2(f) -- lead term from e_2*e_3
```

### Step 2:
Induce valuation on \(\mathbb{C}[x_1,x_2,x_3]^{A_3}\)
```
vA0 = valM(R, v0);
vA1 = valM(R, v1);
vA2 = valM(R, v2);
```

### Example:
```
use R;
g = x_1^2 + x_2^2 + x_3^2
vA0(g)
vA1(g)
vA2(g)
```
```
h = (x_1^2 - x_2^2)*(x_1^2 - x_3^2)*(x_2^2 - x_3^2)
vA0(h)
vA1(h)
vA2(h)
```


## Future Directions
* General development
* Use valuations in other M2 pacakges:

1. [Tropical](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/Tropical/html/index.html) -- tropical computations
2. [FormalGroupLaws](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/FormalGroupLaws/html/index.html) -- commutative formal group laws
3. [SubalgebraBases](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.18/share/doc/Macaulay2/SubalgebraBases/html/index.html) -- canonical subalgebra bases (Sagbi bases)
* Puiseux series
* Suggestions?
