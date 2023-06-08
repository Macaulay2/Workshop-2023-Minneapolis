# Valuations
## Introduction and Motivation
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
   (-) [Tropical](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/Tropical/html/index.html) --  the main M2 package for tropical computations
   (-) [FormalGroupLaws](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.17/share/doc/Macaulay2/FormalGroupLaws/html/index.html) -- commutative formal group laws
   (-) [SubalgebraBases](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.18/share/doc/Macaulay2/SubalgebraBases/html/index.html) -- A package for finding canonical subalgebra bases (Sagbi bases)
3. Puiseux series
4. Any suggestions?