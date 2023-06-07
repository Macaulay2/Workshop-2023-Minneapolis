
-- The function diagonalizeChangeBasisMatrix takes a symmetric non-singular matrix A and 
-- returns (D, Q), where D is a diagonal matrix and Q is the change of basis matrix such that
-- D= QAQ^T

-- The reduction does not use any divisions and works over the overall ring that A
-- is defined.  For example, it will do the reduction over ZZ. 

-- Example:  A = matrix{{1,2},{2,5}}
-- diagonalizeChangeBasisMatrix A returns  (D, Q),  D=((1,01),(0,1)),  Q=((1, 0), (-2,1))

diagonalizeChangeBasisMatrix = method()

diagonalizeChangeBasisMatrix (Matrix) := (Matrix, Matrix) => (A) -> (
-- input is  a nxn matrix (A)   
    R:= ring A;
    n := numRows(A);
    I := id_(R^n);
    m:= matrix A;
   
-- We create the matrix (a|I)
    N := m|I;

-- don't need    M := mutableMatrix N;

-- if n=1, it will be returned.      
    if (n==1) then (
	return (1, N);
    )
    else(
	-- if n>1 we will send it to the diagonalize2 routine to return it in the form
	-- (D | Q), where D = Q A Q^T 
	ans := matrix diagonalizeCB(N);
	D:=submatrix'(ans, , {n..(2*n-1)} );
	Q:=submatrix'(ans, ,{0..(n-1)} );
	return (D, Q);
	);
);

-- Function used internally
-- Input:  M a n x k matrix
-- Output: N a n x k matrix where the first n columns have been row/column reduced to a diagonal matrix
-- The function works by using row/column operations to clear entries in the first row/column (except 
--    for a11) adn then recursively apply program to smaller (n-1) x (k-1) matrix.  

--- Example:  A = matrix{{1,2},{2,5}},  I = 2 x2 identity,  M= (A I). 
-- diagonalizeCB A outputs  (I Q), with Q=((1,0), (-2, 1)).

diagonalizeCB = method()

diagonalizeCB (Matrix) := (Matrix) => (N) -> (

--input: M is a n x 2n matrix.  M needs to be reduced.     
    M := mutableMatrix N;
    n := numRows M;
    if (M_(0,0) == 0) then (	
	--  We first want to see if all diagonal entries are 0.  if so, need to create a non-zero diagonal element
	flag := 1;
	row := 0;
	while flag > 0 and row < n-1 do (
	    row=row +1;
	    if M_(row,row) != 0 then (
		 --we have found non-zero entry
		 flag =0;
		 );
            );
	-- if flag = 0 then there is a non-zero diagonal value
	if flag == 1 then (
	    -- find first entry in 0th row that is non-zero
	    flag2 := 1;
	    row = 0 ;
	    while flag2 == 1 and row < n-1  do (
		row = row +1;
		if  (M_(0, row) ) != 0 then flag2 = 0;
		);
	    -- if flag2 =0, then row gives the row/coln number with A_(0,row) nonzero
	    if flag2 == 0 then  (
		-- we now create a matrix with a non-zero diagonal element
		rowAdd(M, 0, 1, row);
	    	columnAdd(M, 0, 1, row);
		)
	    else (
		print "Error: Matrix A was singular"; 
		return N;
		);
	    );
	);
   
    --  We then analayze the case when we know that one of the diagonal elements is non-zero
    flag = 1;
    row = 0;
    while flag > 0 and row < n-1 do (
	if M_(row,row) != 0 then (
	    flag =0;
	    if row >0 then  (
	    --we have found non-zero entry
	       rowSwap(M,row,0);
	       columnSwap(M,row,0);
	       );
	    );
	 row=row +1;
	);
    if flag == 1 then (
	print "Error: Original square matrix A was singular"; 
	return;
	);
    
    
    --entry in M_(col,col) is non-zero at this point
    for row from (1) to (n-1) do (
	temp:=M_(row,0);
	rowMult(M, row, M_(0,0));
	rowAdd(M, row,-temp,0);	
	);
   
    for col from 1 to (n-1) do (
	temp := M_(0,col);
        columnMult(M, col, M_(0,0));		
	columnAdd(M, col,-temp, 0);	
	);
   
    
    -- At this point, we have taken care of first row, column; 
    -- First part of M is ( a11 0 |  0 A')
    if n==2 then (
	-- if n=2, then we have the final matrix version (D|Q) and can return it
	return matrix M; )
    else(
	-- if n>2, then we need to strip the top row and column from M, then send the new matrix
	-- to this routine
	M1 := matrix M;
	toprow := M1^{0};
	firstcolumn := M1_{1};
        firstcolumnTr := submatrix'(firstcolumn, {0}, );
	-- gets rid of first row and first column
	N1 := submatrix'(M1, {0},{0});
	M1 = diagonalizeCB(N1);
	-- Once we return, we reconstruct the final matrix from M1 and return it.
	M3 := firstcolumnTr| M1;
	M4 := toprow || M3;
	M5 := matrix M4;
	return M5;
	) ;
);


--  Input:matrix (A, Q), 
--  Output QAQ^T
--  Example:  A = matrix{{1,2},{2,5}},  Q= matrix {{1,0},{-2,1}},  
--  checkQFactorization (A, Q) outputs  ((1,0),(0,1))

checkQFactorization = method()

checkQFactorization (Matrix, Matrix) := (Matrix) => (A, Q) -> (
    Q2:= transpose Q;
    return Q*A*Q2;
    );

-- Check that the diagoanlization process worked
-- Input matrix A
-- Output (A, D, Q, Q*A*Q^T), where D diagonal, Q change of basis matrix 
-- If QA(Q^T)<> D, then output error

checkQFactorizationAll = method ()
 
checkQFactorizationAll (Matrix) := (Matrix, Matrix, Matrix, Matrix) => (A) -> (
    B:= diagonalizeChangeBasis(A);
    if B_0 == (B_1)*A*(transpose B_1) then  (	
         return (A, B_0, B_1, (B_1) * A * (transpose B_1),"Diagonalization works");
	 ) else (
	 return  (A, B_0, B_1, (B_1) * A * (transpose B_1), "Diagonalization failed");
    );
);
    
