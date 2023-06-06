load "ABuilder.m2"
load "DegenMap.m2"
faceMapi = (n,k,i) -> (
    Abigvector = apply(ABuilder(n,k), row->row_i);
    Asmallvector = apply(ABuilder(n-1,k), row->row_(i-1));
    outputMat = new MutableList from binomial(n-1,k):(new List from binomial(n,k):0);
    inputMat = listBuilder(binomial(n,k));
    sumCounter = 0;
    sumPositions = positions(Asmallvector, i->i>=1);
    zeroCounter = 0;
    for l to #inputMat-1 do (
    	if Abigvector#l == 0 then (outputMat#zeroCounter = inputMat#l;
	    zeroCounter=zeroCounter+1
	)
	else if Abigvector#l==1 then (
	    myIndex = sumPositions_sumCounter;
	    outputMat#myIndex = outputMat#myIndex + inputMat#l;
	    sumCounter=sumCounter+1;
    	);
    );
    toList outputMat
)
    
