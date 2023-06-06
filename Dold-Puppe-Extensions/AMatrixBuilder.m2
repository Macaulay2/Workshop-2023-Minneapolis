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

ABuilder = (n,k) -> (
    L := sort compositions(k+1,n-k);
    matrix (
    for mu in L list
    	for j to n list
    	    entryCalculator(mu,j)
    )
)	

AhashBuilder = (n,k) -> (
    L := sort compositions(k+1,n-k);
    H := new MutableHashTable;
    for i from 0 to #L-1 do (
	for j to n do (H#(i,j) = entryCalculator(L_i,j)));
    H
    ) 	


