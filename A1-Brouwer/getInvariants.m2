load "./diagonalize.m2"
loadPackage "RationalPoints2"

getInvariants=method()

getInvariants (Matrix) := (List) => (A) -> (
    n:= numRows(A);
    kk:=ring A;

    discrRep:= det A; --representative of dicriminant of bilinear form

    if (kk===RR or instance(kk,RealField) or kk===QQ) then(
        diagA := diagonalize(A);
        posEntries := 0; --for-loop counts the number of positive diagonal entries of diagA
         for i from 0 to (n-1) do(
            if diagA_(i,i)>0 then(
                posEntries=posEntries+1;
            );
        );

        negEntries := n-posEntries;
        signature:=posEntries-negEntries;

        return {n,signature,discrRep}) --output signature if working over ordered fields

    else return {n,discrRep}; -- don't output signature if not working over ordered fields
);



D=matrix{{-1/1,-1,1,1},{-1,1,1,0},{1,1,0,0},{1,0,0,0}}
print diagonalize(D)
print getInvariants(D)

isIdenticalDiscriminant=method()

isIdenticalDiscriminant (Matrix,Matrix) := (Boolean) => (A,B) -> (
    a:=det(A);
    b:=det(B);
    kk:= ring A;
    R:=kk[x];

    f= x^2 - (a/b);

    L:=zeros(f)






);



