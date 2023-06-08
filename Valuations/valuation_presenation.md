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

## Ordered $Q$-modules
## Kaveh-Manon Example
## Future Directions
