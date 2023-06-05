# Alternating Sign Matrix Varieties (ASMs).

Two very important results underly a great deal of computation in this package.  The first is that the Fulton generators of ASM ideals form a Gröbner basis (Knutson-Miller+Weigandt) under all antidiagonal term orders.  Consequently, the (unique) antidiagonal initial ideals of ASM ideals are radical (as, therefore, are ASM themselves). The other key result is that a homogeneous ideal shares its Castlenuovo-Mumford regularity (henceforth, regularity) and depth with all of its radical initial ideals, should it have any (Conca-Varbaro).  Together, these results tell us that, in order to compute the regularity and depth of an ASM ideal, it suffices to compute the same for its antidiagonal initial ideal, which is at least as easy to write down as the ASM ideal itself and typically much faster to resolve. An ASM ideal and its antidiagonal initial ideal will also share the same multidegree (in particular, the same single and double Schubert polynomial).

The basic functions in this package will (listed by goal, not by individual function name) 
⭐ = "priority for workshop"

1. ✔ Take a permutation (as a list or as a matrix) or ASM or rank table and produce the associated ASM ideal (via Fulton generators) or antidiagonal initial ideal (via initial terms of Fulton generators). **[requires tests and documentation]**
2. ✔ Take a permutation (as a list or as a matrix) or ASM and produce its (minimal) rank table. **[requires tests and documentation]**
3. ⭐ [beginner] Take a minimal rank table and produce its ASM.  See Weigandt, Lemma 3.2.
4. ✔ Take a permutation (as a list or as a matrix) or ASM or rank table and produce its Rothe diagram (meaning set of ordered pairs, not visualization). **[requires tests and documentation]**
5. ⭐ [beginner] Give the height of an Schubert determinantal ideal by taking the length of the permutation/counting the size of its diagram. (Caution: this is only guaranteed to give the correct height for Schuberts, not all ASMs.)
6. ⭐ [beginner] Give the height, depth, and regularity of an arbitrary ASM ideal by taking the same of its antidiagonal initial ideal.
7. ✔ Take a permutation (as a list or as a matrix) or ASM or rank table and produce its essential set. **[requires tests and documentation]**
8. ⭐ [beginner] Take a permutation (as a list or as a matrix) and produce its single or double Schubert polynomial and single or double Grothendieck polynomial using divided differences.
9. Convert a permutation among its list, permutation matrix, and bijection forms.

More complicated general purpose functions:

11. ⭐ Take an arbitrary ASM and produce its multidegree for the usual N^n or N^{2n} grading via its antidiagonal initial ideal using degeneration and additivity of multidegrees.
12. ⭐ Take an ASM ideal and find its prime decomposition via its antidiagonal initial ideal: Each component of the antidiagonal initial ideal indexes a (not-necessarily-distinct) Schubert determinantal ideal in the decomposition of the ASM ideal itself.  The generators each such component are variables whose indices give a reduced word for one such permutation (Knutson-Miller).

All ASM varieties are unions of matrix Schubert varieties; however, not every union of matrix Schubert variety is an ASM variety.  For this reason, it comes up with some frequency that we have an ideal in hand and want to know if it is an ASM ideal.

13. ⭐ Decide whether or not I is an intersection of Schubert determinantal ideals: Take the antidiagonal initial ideal of I and verify that it is radical.  Assuming it is, read the set of permutations, Perm(I), from the components of the antidiagonal initial ideal of I.  (That is, apply the function written for 12.)  Ask whether the intersection of Schuberts indexed by Perm(I) is indeed I.

13.5. Suppose we know by 13 that I is an intersection of Schubert determinantal ideals.  Take mins of rank tables of the Schuberts involved, build the ASM ideal J of that rank table and then check if I=J.

The set of ASM ideals is closed under addition.  Typically it's not easy to see which ASM is the matrix corresponding to a given sum of ASM ideals.

14. ⭐ Given a set of ASMs, find the rank table for their sum by taking entrywise minima of the rank table of each of the given ASMs:  Find the ASM of the output rank table (using the function written for 3, rank table will be minimal).
15. Given a non-minimal rank table, find its ASM: Exclude the vacuous rank conditions and assign the appropriate bigrassmannian to the others.  Apply the function from 14 on those bigrassmannians.

Regularity:

16. Compute regularity of any CM ASM (allow opt out of CM check, opt our automatically if given a permutation) as degree of K-polynomial minus height
17. Decide if a permutation is vexillary.
18. Decide if a permutation is 1432 avoiding.
19. Implement Theorems 1.1 and 1.5 of Rajchgot-Robichaux-Weigandt.
20. ⭐ Compute Rajchgot index of a permutation and implement Theorem 1.1 of Pechenik-Speyer-Weigandt.

While we're at it, if time: 

21. Decide if a permutation is CDG (mild generalization of function from 16).
22. Give the (unique) diagonal initial ideal of a CDG permutation.  See Klein, Conjecture 1.1 and paragraph before it.
23. Decide if a permutation is Cartwright-Sturmfels (mild generalization of function from 16) and, if so, give its universal Gröbner basis.  See Conca-De Negri-Gorla Theorem 5.6.
