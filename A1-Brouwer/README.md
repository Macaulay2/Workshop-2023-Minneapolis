# A1-Brouwer degrees

The goal of this project is to implement algorithms for computing $\mathbb{A}^1$-Brouwer degrees, valued in the Grothendieck-Witt ring of symmetric bilinear forms over a field. We'll be following work in [*Bézoutians and the A1-degree* by B., McKean, Pauli](https://www2.math.upenn.edu/~tbraz/bezoutian.pdf), abbreviated [BMP].

Our goals for the workshop:

1. Implement the Grothendieck-Witt ring into Macaulay2. For background see [GW.pdf](./References/GW.pdf)
2. Implement and improve upon the algorithms from [BMP] for computing $\mathbb{A}^1$-degrees. For background see [Bezoutians.pdf](./References/Bezoutians.pdf)

The other links may also be useful:
- [Diagonalizing symmetric bilinear forms](https://www2.math.upenn.edu/~tbraz/notes/diagonalizing-forms.pdf), we could potentially use this to exploit sparsity in Gram matrices to more quickly diagonalize forms. The Gram matrices attached to A1-degrees often have this sort of structure, so it would be highly useful.
- [Localization as quotients](https://machineappreciation.wordpress.com/2022/04/15/localizations-as-quotients/), a blog post about normal bases and saturation.
