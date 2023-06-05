entryCalculator = (mu,j) -> (
    u=0;
    for i from 0 to #mu-1 do (
	u=u+mu_i;
	if u == j+1 then (
		if mu_i == 1 then
			return 2
		else 
			return 1
		);
	);
    0
)
	
