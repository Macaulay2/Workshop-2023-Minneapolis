# A1-Brouwer degrees

The goal of this project is to implement algorithms for computing $\mathbb{A}^1$-Brouwer degrees, valued in the Grothendieck-Witt ring of symmetric bilinear forms over a field. We'll be following work in [*Bézoutians and the A1-degree* by B., McKean, Pauli](https://www2.math.upenn.edu/~tbraz/bezoutian.pdf), abbreviated [BMP].

Our goals for the workshop:

1. Implement the Grothendieck-Witt ring into Macaulay2. For background see [GW.pdf](./References/GW.pdf)
2. Implement and improve upon the algorithms from [BMP] for computing $\mathbb{A}^1$-degrees. For background see [Bezoutians.pdf](./References/Bezoutians.pdf)

The other links may also be useful:
- [Diagonalizing symmetric bilinear forms](https://www2.math.upenn.edu/~tbraz/notes/diagonalizing-forms.pdf), we could potentially use this to exploit sparsity in Gram matrices to more quickly diagonalize forms. The Gram matrices attached to A1-degrees often have this sort of structure, so it would be highly useful.
- [Localization as quotients](https://machineappreciation.wordpress.com/2022/04/15/localizations-as-quotients/), a blog post about normal bases and saturation.


# Completed
✅ Made `GrothendieckWittClass` type with a constructor out of a matrix over a field, with `isWellDefined` boolean asking whether it defines a symmetric bilinear form over a field of characteristic not equal to 2.

✅ Made `gwAdd` and `gwMultiply` methods for adding and multiplying classes in GW

✅ Made `diagonalize` method for diagonalizing symmetric bilinear forms

✅ Made `rationalSimplify`, which takes in a symmetric bilinear form over $\mathbb{Q}$ and diagonalizes it, replaces each diagonal entry by the smallest magnitude integer in its square class, and splits off occurences of $\langle a\rangle + \langle -a\rangle$ as hyperbolic forms

✅ Made `easyIsomorphicGW` which checks, up to a certain height, whether two Gram matrices are congruent.

✅ Added `isIsomorphic2` to test whether two Grothendieck-Witt classes are isomorphic; only implemented over $\mathbb{R}$ and $\mathbb{C}$

# In progress

🔸 Implement a "strategy" for the `wittDecomp` algorithm for detecting the Witt class of a symmetric bilinear form

🔸 Refined algorithms for quickly diagonalizing sparse matrices (see `easyUpperLeftTriangular`)

🔸 Check whether some of our matrix algorithms extend for singular cases (for example `isIsomorphic2`).

🔸 Add tests using `assert(-)` statement

🔸 Documentation, written in M2 language


# Still to do

🔺 Implement A1-Brouwer degrees into Macaulay2

🔺 Figure out how to send the diagonal form of a GrothendieckWittClass to the `CacheTable` in the type file
