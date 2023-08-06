installPackage("A1BrouwerDegrees")
viewHelp A1BrouwerDegrees


M = matrix(QQ,{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
beta = gwClass(M)
integralDiagonalRep(beta)

installPackage("A1BrouwerDegrees", RerunExamples=>true)
check A1BrouwerDegrees


T4 = CC[v];
f5 = {v^4 + v^3 - v^2 - v};
f5GD = globalA1Degree(f5);
f5GD




T4Q = QQ[v];
f5Q = {v^4 + v^3 - v^2 - v};
f5QGD = globalA1Degree(f5Q);
f5QGD


ff = GF(17);
T3 = ff[y_1..y_3];
f3 = {y_1^2, y_2^2, y_3^2};
f4 = {y_2^2, y_3^2, y_1^2};
f3GD = globalA1Degree(f3)
f4GD = globalA1Degree(f4)


ff = GF(17);
T3 = ff[w_1..w_3];
f3 = {w_1^2, w_2^2, w_3^2};
globalA1Degree(f3)


M = matrix(GF(13),{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
beta = gwClass(M);
c=simplifyForm(beta)

beta
M
beta.cache.diagonalForm

class c



-- Check if a matrix is square
isSquare = method()
isSquare (Matrix) := Boolean => M -> (
    numRows(M) == numColumns(M)
)

-- Check if a matrix is diagonal
isDiagonal = method()
isDiagonal (Matrix) := Boolean => M -> (
    if not isSquare(M) then error "Error: matrix isn't square";
    n := numRows(M);
    for i from 0 to n-1 do(
	for j from 0 to n-1 do(
    	    if i != j then(
		if not M_(i,j) == 0 then(
		    return false
		    );
		);
	    );
	);
    true
    )


N = matrix(QQ,{{1,0,0},{0,1,0},{0,0,1}})
N1 = matrix(QQ,{{1,1,0},{0,1,0},{0,0,1}})
N2 = matrix(QQ,{{1,0,1},{0,1,0},{0,0,1}})
N3 = matrix(QQ,{{1,0,0},{1,1,0},{0,0,1}})
N4 = matrix(QQ,{{1,0,0},{0,1,1},{0,0,1}})
N5 = matrix(QQ,{{1,0,0},{0,1,0},{1,0,1}})
N6 = matrix(QQ,{{1,0,0},{0,1,0},{0,1,1}})

isDiagonal(N)
diagonalize(N)
