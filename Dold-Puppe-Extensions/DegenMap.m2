load "ABuilder.m2"
DegenMapHelper = (n, k, i) -> (
    curList = entries id_(ZZ^(binomial(n-1,k)));
    zeroes = new List from (binomial(n-1,k)) : 0;
    A = ABuilder(n,k);
    sig = {};
    l = 0;
    for j to #A-1 do(
	if (A_i_j == 0) then (
	    sig = append(sig, curList_l);
	    l = l+1;
	    continue;
	    ); 
	sig = append(sig, zeroes);
	);
    sigout = matrix(sig);
    sigout
    )
