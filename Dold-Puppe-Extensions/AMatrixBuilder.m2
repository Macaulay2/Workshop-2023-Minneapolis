load "Compositions.m2"

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
	
