Additional code that should be written:
1. A proper cyclotomic polynomial code --started on
2. The nullhomotopy commands are implicitly assuming
   the factorization has length 2 (the definition of a nullhomotopy
   for a longer factorization is different)
3. It would be nice if the branched cover code actually worked
4. Homology of a ZZ/2-graded factorization --DONE
5. The collapse command for longer factorizations should be working properly
6. An "isZZdComplex" command which just checks if the differentials
   all compose to 0 --DONE
7. hh command and eulerChi command (Length of homology and Euler characteristic)
8. An "Unfold" command that converts a factorization into a complex 
   (this makes me think of the test: Unfold a factorization of f, -- DONE
       tensor with S/(f), then check that this is a well-defined complex)
9. Code for the suspension of a factorization
10. Code "trivialFactorization" that takes as inputs a module, an integer,
    and a ring element and outputs the "trivial" factorization (f,1,1 ... ,1)
11. Make function randomTailMF which takes in a polynomial (and upper bound 
    for presentation matrix) and generates a random matrix (module pres) over hypersurface 
    and applies tailMF
