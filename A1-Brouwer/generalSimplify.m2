path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");
needs "GW-type.m2"
load "diagonalize.m2"
load "diagonalForm.m2"
load "diagonalizeOverInt.m2"
load "easyIsomorphicGW.m2"
load "isIsomorphic2.m2"
load "getInvariants.m2"
load "radical.m2"
load "rationalsimplify.m2"
load "splitoffobvioushyperbolics.m2"
load "simplifyForm.m2"
load "squarefreepart.m2"
load "wittDecomp.m2"
load "wittDecompGeneral.m2"

gwSimplify=method();

--This method returns an isomorphic GW class with a nice diagonal representative.
gwSimplify (GrothendieckWittClass) := (GrothendieckWittClass, ZZ, ZZ) => (alpha) -> (
    k=baseField(alpha);
    if (k === RR or instance(k,RealField)) then (
        (wittIndex, anisotropic) := wittDecompInexact(alpha.matrix);
        return (diagonalForm(alpha), wittIndex, numRows(anisotropic));
    );
--When the base field is RR, diagonalForm gives the nicest possible diagonal form
--and wittDecompInexact gives the Witt Index. Dimension of anisotropic part is given
--by number of rows of anisotropic part, also given by wittDecompInexact.

    if (k === CC or instance(k,ComplexField)) then (
     return (diagonalForm(alpha)), 0, numRows(alpha.matrix);
     );

--When the base field is RR, diagonalForm gives the form represented by the identity 
--matrix, with wittIndex 0.

    if (k===QQ) then (
        output=mutableMatrix(QQ, {{}});
        A=diagonalize(alpha.matrix);
        B = rationalSimplify(A);
        --The anisotropic part after factoring out the squares
        (wittdim, mat):= wittDecomp(alpha.matrix);
        n=0;
        for i from 0 to (wittdim-1) do (
            output=safeBlockSum(output, matrix(QQ, {1/1, 0}{0, -1/1}));
        );
        output=safeBlockSum(output, B);
        C = matrix(output);
        return (gwClass(C), wittdim, numRows(B));
    );
    if (instance (k, GaloisField)) then (
        (resultGF, str) := simplifyForm(alpha.matrix);
        return (resultGF, wittIndexGF(alpha), (numRows(alpha.matrix)-2*num));
    );
);

--Witt Index for Finite Field case (extracted from simplifyForm.m2):
wittIndexGF=method()
wittIndexGF (GrothendieckWittClass) := ZZ => (beta) -> (
    if (instance (baseField(beta), GaloisField)===false) then (print "error: base field of input is not a finite field";)
    else (
        diag := {};
	for i from 0 to sub((n-1),ZZ) do(
	    diag = append(diag, A_(i,i));
	    );
        Squares = 0;
	    NonSquares = 0;
	nonSquareRep := sub(0,k);

	for x in diag do(
	    if legendreBoolean(x) then(
		Squares = Squares + 1;
		);

	    if not legendreBoolean(x) then(
		nonSquareRep = x;
		NonSquares = NonSquares + 1;
		);
	    );
	
	
	if legendreBoolean(sub(-1,k)) then(
            wittGFSquare := floor(Squares/2) + floor(NonSquares/2);
            return wittGFSquare;
    );
	
	   if not legendreBoolean(sub(-1,k)) then (
	       wittGFNonSquare := min(numSquares,numNonSquares);
           return wittGFNonSquare;
       ); 
    );
);


--This method returns a diagonal matrix that is congruent to the input.
matrixSimplify=method(
    Options => {
	HeightBound => 4  
	}
)
--matrixSimplify (matrix) := (matrix) => opts -> (gamma) -> (
    --Still writing this part
--);