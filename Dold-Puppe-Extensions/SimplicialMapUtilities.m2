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

RowOfId = (n,l) -> splice {l:0,1,(n-1-l):0}    

FaceMap = (n,k,i) -> (
    Abigvector = apply(ABuilder(n,k), row->row_i);
    Asmallvector = apply(ABuilder(n-1,k), row->row_(i-1));
    myIndices = 0..(binomial(n,k)-1);
    filter = partition(myInd -> Abigvector#myInd,myIndices);
    outputMat = new MutableList from filter#0 / (l -> rowOfId(binomial(n,k),l));
    sumPositions = positions(Asmallvector, i->i>=1);
    for l to #sumPositions-1 do (
	myIndex = sumPositions_l;
	outputMat#myIndex = outputMat#myIndex + (rowOfId(binomial(n,k),(filter#1)#l));
    );
    toList outputMat
)

   
DegenMap = (n, k, i) -> (
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
