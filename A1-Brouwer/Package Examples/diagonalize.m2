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
            if A_(col,col)==0 then (print "Error: Matrix A was singular"; return A;);
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
--wittDecompFromList method & wiitDecompFromMatrix method
wittDecompFromList =method()
wittDecompFromMatrix =method()

wittDecompFromMatrix (Matrix,Ring) := (ZZ,List) => (A,k) -> (
    D:= diagonalize(A);
    return wittDecompFromList((for i from 0 to numRows(A)-1 list (D_(i,i))),k);
)

wittDecompFromList (List,Ring) := (ZZ,List) => (vals,k) -> (
    n:=#vals;
    R:=k[x_0..x_(n-1)];
    f:=sum (for i from 0 to n-1 list (vals_i*x_i^2));
    use k;
    soln:=new MutableList;
    solnPt:=new MutableList;
    for bound from 1 to 10 do (
        soln:= new MutableList from rationalPoints(ideal(f),Bound=>bound);
        if (#soln >1) then(
            if soln#0 == toList(n:0) then (solnPt=soln#1) else (solnPt=soln#0);
            break;
        );
    );
    y:= new MutableList from toList(n:0);
    if solnPt#0==0 
        then (y#0=1;) 
        else (y#0=-solnPt#1; y#1=solnPt#0);
    --solnPt and y span a copy of |H in A
    --now turn solnPt and y into row matrices
    v1:=matrix{toList solnPt};
    v2:=matrix{toList y};
    print v2;
    return 0;
);

print wittDecompFromMatrix(B,QQ);



