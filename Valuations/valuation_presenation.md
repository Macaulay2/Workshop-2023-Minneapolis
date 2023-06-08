# Valuations
## Setup and Motivation

### What is a valuation?

A valuation is a map \(v: R \rightarrow \Gamma\) between a ring \(R\) and an ordered semigroup \(\Gamma\) such that:

--* \(v(ab) = v(a) + v(b)\)
--* \(v(a+b) \ge \min\{v(a), v(b)\}\)
--* (For \(k\)-algebras) \(v(\lambda a) = v(a)\)

They are fundamental for the definition of:

--* Tropical varieties (polyhedral complexes that encode data about the original variety)
--* Khovanskii bases (a generating set for an algebra analogous to Groebner bases)

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
>>>>>>> dbf396aef443d8af0de47799471ab842a2b4b2df
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
Code
```

### Step 2:
Induce valuation on \(\mathbb{C}[x_1,x_2,x_3]^{A_3}\)
```
Code
```

### Example:
```
Code
```


## Future Directions
1. General development
2. Use valuations in other M2 pacakges:
   - [Tropical](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/Tropical/html/index.html) --  the main M2 package for tropical computations
   - [FormalGroupLaws](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/FormalGroupLaws/html/index.html) -- commutative formal group laws
   - [SubalgebraBases](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.18/share/doc/Macaulay2/SubalgebraBases/html/index.html) -- A package for finding canonical subalgebra bases (Sagbi bases)
3. Puiseux series
4. Any suggestions?
