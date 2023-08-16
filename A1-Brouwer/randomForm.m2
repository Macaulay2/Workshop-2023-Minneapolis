needs ("GW-type.m2")
--This file generates a random GW-class represented by a matrix 
--of given rank n over a given field k

randomForm=method()

randomForm(InexactFieldFamily, ZZ) := (GrothendieckWittClass) => (K, n) -> (
    outputMatrix=mutableMatrix(id_(K^n));
    B=0;
while B===0 do (
    for i from 0 to n-1 do (
        for j from i to n-1 do(
            if (K===QQ) then (randomNum=random QQ_100)
            else if (K === RR or instance(K,RealField)) then (randomNum=random RR_100)
            else if (K === CC or instance(K,ComplexField)) then (randomNum=random CC_100);
            outputMatrix_(i, j)=randomNum;
            outputMatrix_(j, i)=randomNum;
        )
    );
    B=det outputMatrix;
);
    return gwClass(matrix(outputMatrix));
);

randomForm(GaloisField, ZZ) := (GrothendieckWittClass) => (K, n) -> (
    outputMatrix=mutableMatrix(id_(K^n));
    B=0;
while B===0 do (
    for i from 0 to n-1 do (
        for j from i to n-1 do(
            randomNum=random K;
            outputMatrix_(i, j)=randomNum;
            outputMatrix_(j, i)=randomNum;
        )
    );
    B=det outputMatrix;
);
    return gwClass(matrix(outputMatrix));
);