needsPackage "RationalPoints2"
--diagonalize method
diagonalize = method()

diagonalize (Matrix) := (Matrix) => (AnonMut) -> (
    A := mutableMatrix AnonMut;
    n:=numRows(A);
    for col from 0 to (n-1) do (
        if A_(col,col) == 0 then (
            for row from col+1 to n-1 do (
                if A_(row,col) != 0 then (
                    --we have found non-zero entry
                    rowAdd(A,col,1,row);
                    columnAdd(A,col,1,row);
                    break;
                );
            );
            if A_(col,col)==0 then (print "Error: Matrix was singular"; return A;);
        );
        --entry in A_(col,col) is non-zero at this point
         for row from (col+1) to (n-1) do (
            temp:=A_(row,col);
            rowAdd(A,row,-temp/A_(col,col),col);
            columnAdd(A,row,-temp/A_(col,col),col);
        );

    );
    return matrix A
)

A=matrix{{1/1, 0},{0 ,1}};
B=matrix{{0,1/1},{1 ,0 }};
C=matrix{{1,1/1},{1 ,1 }};
D=matrix{{-1/1,-1,1,1},{-1,1,1,0},{1,1,0,0},{1,0,0,0}};

print diagonalize(A)
print diagonalize(B)
print diagonalize(C)
print diagonalize(D)

print "------------------------------------------------------------";
--wittDecomp method

wittDecomp =method()
wittDecomp (Matrix,Ring) := (ZZ,Matrix) => (A,k) -> (
    n:=numRows(A);
    R:=k[x_0..x_(n-1)];
    f:=sum (
        for i from 0 to n-1 list (
            sum (for j from 0 to n-1 list (A_(i,j)*x_i*x_j))
            )
        );
    use k;
    solnPt:=new MutableList;
    solnFound = false;
    for bound from 1 to 10 do (
        solns:= new MutableList from rationalPoints(ideal(f),Bound=>bound);
        if (#solns >1) then(
            if solns#0 == toList(n:0) then (solnPt=solns#1) else (solnPt=solns#0);
            solnFound =true;
            break;
        );
    );
    --if no solutions found we assume the form is anisotropic
    if ( not solnFound) then (return (0,A));
    --if solution found for rank 2 form, then the form is purely isotropic
    if (n==2) then (return (1,{}));

    --find y not orthogonal (wrt bilinear form) to x
    x:=matrix{toList solnPt}; --x as a row matrix
    xA:=x*A; --x*A
    y :=new MutableMatrix from matrix{toList(n:(0/1))};
    for i from 0 to (n-1) do (
        if (xA_(0,i) != 0) then (y_(0,i)=1; break;);
    );
    --now x and y span a copy of |H in the bilinear form
    --we need to find a basis of vectors orthogonal (wrt bilinear form) to x and y
    Red := reducedRowEchelonForm (x||matrix(y));
    W := mutableMatrix matrix(toList((n-2):toList(n:0/1))); --W will contain as its rows w2,...,w(n-1) orthogonal to x and y
    for i from 2 to (n-1) do (
        W_(i-2,i)=1;
        W_(i-2,0)=-Red_(0,i);
        W_(i-2,1)=-Red_(1,i);
    );
    --now recursively apply wittDecomp to W*A*W^T a (n-2)-by-(n-2) Gram matrix
    Wmat := matrix(W);
    subComputation := wittDecomp(Wmat*A*transpose(Wmat),k);
    return (1+subComputation_0, subComputation_1);
)

print wittDecomp(D,QQ);


