
--load ("/home/macaulay/A1-Brouwer/diagonalize.m2")
--load ("/home/macaulay/A1-Brouwer/easyIsomorphicGW.m2")
load "diagonalize.m2"
load "GW-Type.m2"
load "easyIsomorphicGW.m2"

diagonalize = method()

diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    k := ring AnonMut;
    if isField k == false then error "Error: expected matrix entries from a field";
    A := mutableMatrix AnonMut;
    if A != transpose(A) then (
        error "Matrix is not symmetric";
	);
    n := numRows(A);
    for col from 0 to (n-1) do (
	--If diagonal entry in column "col" is zero
        if A_(col,col) == 0 then (
            for row from col+1 to n-1 do ( 
		--Scan for nonzero entries below the diagonal entry
                if A_(row,col) != 0 then (
		    --Row reduction to make A_(col,col) nonzero
                    rowAdd(A,col,1,row);
		    --Column reduction to keep reduced matrix congruent to original matrix
                    columnAdd(A,col,1,row);
                    break;
                );
            );
        );
        --If nonzero entry at or below A_(col,col) is found, we use it to clear the column below
        if (A_(col,col)!=0) then (
            for row from (col+1) to (n-1) do (
                temp:=A_(row,col);
		--More row reduction make every entry below A_(col,col) is zero
                rowAdd(A,row,-temp/A_(col,col),col);
		--Column reduction to keep reduced matrix congruent
                columnAdd(A,row,-temp/A_(col,col),col);
            );
        );

    );
    return matrix A 
    )

GrothendieckWittClass = new Type of HashTable
GrothendieckWittClass.synonym = "Grothendieck Witt Class"
matrix GrothendieckWittClass := Matrix => beta -> beta.matrix
gwClass = method()

-- A class in GW can be constructed from a representing matrix
gwClass (Matrix) := GrothendieckWittClass => M -> (
  new GrothendieckWittClass from {
      symbol matrix => M,
      symbol cache => new CacheTable
      }
)
baseField = method()
baseField GrothendieckWittClass := Ring => beta -> (
    ring beta.matrix;
)

easyIsomorphicGW = method(
    Options => {
	HeightBound => 3
	}
    )

easyIsomorphicGW (GrothendieckWittClass, GrothendieckWittClass) := Boolean => opts -> (beta, gamma) -> (
    -- Returns error if the matrices are defined over different base fields
    if not baseField(beta) === baseField(gamma) then error "Error: these classes have different underlying fields";
    
    A := beta.matrix;
    B := gamma.matrix;
    
    -- Return false if the matrices are not of the same size
    if not numRows(A) == numRows(B) then return false;
    
    k:= baseField(beta);
    n := numRows(A);
    
    -- Build a generic matrix P in indeterminants over our field
    R := k[x_(1,1)..x_(n,n)];
    P := genericMatrix(R,n,n);
    
    -- Take the n^2 equations produced by this matrix equality
    testZeroMatx := transpose(P)*A*P - B;
    fullEqns := flatten entries testZeroMatx;
    I := ideal(fullEqns);
    
    -- Increasing the height bound, check for solutions at each stage
    for i in 1..opts.HeightBound do(
	L := rationalPoints(k, I, Bound=>i);
	if not #L === 0 then(
	    return true;
	    break;
	    );
	if #L === 0 then(
	    print("No solutions up to height " | toString(i));
	    );
     	);
    print("Height bound exceeded");
    return false;
    )

testMatrix1 = matrix(QQ,{{6/1,-1/1},{2/1,3/1}});
testMatrix2 = matrix(QQ, {{0,1/1},{1/1, 0}});
resultMatrix1 = matrix(QQ, {{1/1,0},{0,-1/1}});
resultMatrix2 = matrix(QQ, {{-1/1,0},{0,1/1}});
A=diagonalize(testMatrix2);
E=gwClass(A);
B=gwClass(resultMatrix1);
C=gwClass(resultMatrix2);
assert(easyIsomorphicGW(E, B)===true or easyIsomorphicGW(E, C)===true);