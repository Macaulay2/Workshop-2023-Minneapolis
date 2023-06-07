load "wittDecomp.m2"
load "diagonalize.m2"

--wittDecompGeneral returns dimention of totally anistropic part,
--witt index and anisotropic part
wittDecompGeneral = method();

wittDecompGeneral (Matrix) := (ZZ,ZZ,Matrix) => (A)->(
    n=numRows(A);
    k:=ring A;
    D:=diagonalize(A);
    nonZeroDVals := new MutableList;
    for i from 0 to n-1 do (if (D_(i,i)!=0) then (append(nonZeroDVals,D_(i,i))));
    dimReg:=length(toList nonZeroDVals) ;
    regPart := mutableMatrix matrix(toList(dimReg:toList(dimReg:0/1)));
    for i from 0 to dimReg-1 do (regPart_(i,i)=nonZeroDVals_i);
    regWittDecomp := wittDecomp(regPart);
    return (n-dimReg, regWittDecomp_0,regWittDecomp_1);
)

A=matrix{{0/1,0},{0,1}};
print wittDecompGeneral(A);