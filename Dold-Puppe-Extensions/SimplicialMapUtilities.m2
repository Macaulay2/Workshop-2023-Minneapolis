compHash = (n,k) -> (new HashTable from for i in (P=compositions(k,n-k)) list i =>position(P,j->j==i))

entryCalculator = (mu,j) -> (
    u=0;
    for i from 0 to #mu-1 do (
	u=u+mu_i+1;
	if u == j+1 then (
		if mu_i == 0 then
			return 2
		else 
			return 1
		);
	);
    0
)

ABuilder = (n,k) -> table(sort compositions(k+1,n-k),0..n,entryCalculator)

rowOfId = (n,l) -> splice {l:0,1,(n-1-l):0}    

-* Not working yet
faceMapCase0 = (n,len) -> (
    -- If bounded => true then do this
    -- TODO implement option toggle for bounded complexes
    topDeg = sum(len+1,i->binomial(n,i));
    offset = 0 -- keeps track of where we are
    for k to len list (
    -- First loop over the vertical strips of columns
    -- Create appropriately sized 
	for row to topDeg-1 list (
	-- for fixed k, loop over the rows of the matrix
    	    	effectiveRow = row-offset;
		if effectiveRow<binomial(n,k) then splice {row:0,2,(binomial(n,k)-1-row):0}
    	    	-- Put in symbol for differential
		else if (condition) then splice {row:0,1, (binomial(n,k)-1-row):0}
    	    	-- Put in symbol for identity
		else 
		-- put in zeros
    )
*-

faceMapCasen = (n,k) -> (
    Abigvector = apply(ABuilder(n,k), row->row_n);
    Asmallvector = apply(ABuilder(n-1,k), row->row_(n-1));
    onePos = positions(Abigvector, i -> i==1);
    outputMat = onePos / (l -> rowOfId(binomial(n,k),l))
    )
 
faceMapCase1 = (n,k,i) -> (
    Abigvector = apply(ABuilder(n,k), row->row_i);
    Asmallvector = apply(ABuilder(n-1,k), row->row_(i-1));
    zeroPos = positions(Abigvector, i -> i==0);
    onePos = positions(Abigvector, i -> i==1);
    outputMat = new MutableList from zeroPos / (l -> rowOfId(binomial(n,k),l));
    sumPositions = positions(Asmallvector, i->i>=1);
    for l to #sumPositions-1 do (
	myIndex = sumPositions_l;
	outputMat#myIndex = outputMat#myIndex + (rowOfId(binomial(n,k),onePos#l));
    );
    toList outputMat
)

   
degenMap = (n, k, i) -> (
    idSize = binomial(n-1,k);
    zeroes = new List from (binomial(n-1,k)) : 0;
    A = ABuilder(n,k);
    sig = new MutableList from 0..(#A-1);
    l = 0;
    for j to #A-1 do(
	if (A_j_i == 0) then (
	    sig#j = rowOfId(idSize,l);
	    l = l+1;
	    continue;
	    ); 
	sig#j = zeroes;
	);
    new List from sig
    )
