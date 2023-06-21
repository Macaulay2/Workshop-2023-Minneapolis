path = append(path, "/home/macaulay/A1-Brouwer/");
path = append(path, "../A1-Brouwer/");

load "wittDecomp.m2"
load "diagonalize.m2"

--Nikita

--wittDecompGeneral returns dimention of totally anistropic part,
--witt index and anisotropic part
wittDecompGeneral = method();

wittDecompGeneral (Matrix) := (ZZ,ZZ,Matrix) => (A)->(
    n=numRows(A);
    k:=ring A;
    D:=diagonalize(A);
    nonZeroDvals:=new MutableList from {};
    for i from 0 to n-1 do (
        if (D_(i,i)!=0) then (nonZeroDvals = append(nonZeroDvals,D_(i,i));)
    );
    dimReg:=#nonZeroDvals;
    if dimReg!=0 then (
        regPart := mutableMatrix matrix(toList(dimReg:toList(dimReg:0/1)));
        for i from 0 to dimReg-1 do (regPart_(i,i)=nonZeroDvals#i);
        regWittDecomp := wittDecomp(matrix regPart);
        return (n-dimReg, regWittDecomp_0,regWittDecomp_1);
    );
    return (n,0,matrix{{}});
)
