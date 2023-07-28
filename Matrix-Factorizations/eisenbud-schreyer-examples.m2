
needs "ZZdFactorizations.m2"

quadraticMF = method()

--Input: A pair of lists {x_0..x_n}, {y_0..y_n} of elementrs from the same ring
--Output: A matrix factorization for q = \sum x_i*y_i, defined recursively as in Eisenbud--Schreyer

quadraticMF(List, List) := ZZdFactorization => (L, L') -> (
    --First, do some checks
    if not(#L == #L') then error "expected lists of same length";
    S := ring L#0;
    n := #L - 1;
    for i to n do if not instance(L#i, S) or not instance(L'#i, S) then error "expected all list entries to be in the same ring";
    --initialize the recursion
    phi0 := matrix{{L#0}};
    psi0 := matrix{{L'#0}};
    --recursion
    for i from 1 to n do(
	phi = matrix{
	    {L#i * id_(S^(2^(i-1))), phi0},
	    {psi0, (-L'#i)*id_(S^(2^(i-1)))}};
	psi = matrix{
	    {L'#i * id_(S^(2^(i-1))), phi0},
	    {psi0, (-L#i)*id_(S^(2^(i-1)))}};
	phi0 = phi;
	psi0 = psi;
	);
    ZZdfactorization{phi0, psi0}
    )

quadraticMF(Sequence, Sequence) := ZZdFactorization => (X, X') -> (
    quadraticMF(toList X, toList X')
    )

--Helper function
--In: Two lists of ring elements, plus a skew symmetric matrix
--Out: Matrix factorization for quadratic form q as in Eisenbud--Schreyer
quadraticMF(Matrix, List, List) := ZZdFactorization => (Lambda, L, L') -> (
    --Some checks
    if not(#L == #L') then error "expected lists of ring elements of same length";
    n := #L-1;
    if not(rank source Lambda == 2*(n+1)) then error "expected input matrix to have dimension 2N x 2N, where N is length of each input list";
    if not(Lambda + transpose Lambda == 0) then error "expected skew symmetric matrix";
    k := ring Lambda; --should be scalars
    G := matrix{ {0, id_(k^(n+1))}, {id_(k^(n+1)),0}}*Lambda;
    v := flatten entries(matrix{ L|L'}*G);
    quadraticMF(for i to n list v#i, for i to n list v#(n+1+i))
    )

quadraticMF(Matrix, Sequence, Sequence) := ZZdFactorization => (Lambda, X, X') -> (
    quadraticMF(Lambda, toList X, toList X')
    )


--for constructing the Ulrich modules
ulrichFromMF = method()

ulrichFromMF(Matrix, List, List) := Module => (Lambda, L, L') -> (
    F := quadraticMF(L, L');
    F' := quadraticMF(Lambda, L, L');
    S := ring F;
    pres := F.dd_0|F'.dd_0;
    q1 := polynomial F;
    q2 := polynomial F';
    (coker pres)**(S/ideal(q1, q2))
    )

ulrichFromMF(Matrix, Sequence, Sequence) := Module => (Lambda, X, X') -> (
    ulrichFromMF(Lambda, toList X, toList X')
    )

---
--This function should be in the main type file, but it wasn't working for some reason
polynomial = method()
polynomial(ZZdFactorization) := ZZdFactorization => F -> (
    p := F.period;
    comp := product(for i to p-1 list F.dd_i);
    comp_(0,0)
    )


end--
--------------
--Examples
--Note: Error in ZZdFactorizations.m2 means you should just run the functions above before using, rather than loading the file

needs "eisenbud-schreyer-examples.m2"

S = QQ[x_0..y_2]

Q = quadraticMF(x_0..x_2, y_0..y_2)
polynomial Q

Lambda = matrix{{0,1,2,3,4,5}, {-1,0,1,2,3,4},{-2,-1,0,1,2,3},{-3,-2,-1,0,1,2},{-4,-3,-2,-1,0,1},{-5,-4,-3,-2,-1,0}}
Q' = quadraticMF(Lambda, {x_0, x_1, x_2}, {y_0, y_1, y_2})
polynomial Q'

M = ulrichFromMF(Lambda, x_0..x_2, y_0..y_2)
R = ring M
prune M

--Using a skew symmetric matrix coming from diagonal with distinct entries should give ulrich module
use S
D = {4,1,2}
mat = matrix{{0, diagonalMatrix(D)}, {-diagonalMatrix(D), 0}}
N = ulrichFromMF(mat, x_0..x_2, y_0..y_2)
