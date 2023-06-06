load "ABuilder.m2"
listBuilder = n -> for i to n-1 list
    (new List from i:0) | {1} | (new List from (n-1-i):0)

DegenMapHelper = (n, k, i) -> (
    curList = listBuilder(binomial(n-1,k));
    zeroes = new List from (binomial(n-1,k)) : 0;
    A = ABuilder(n,k);
    sig = {};
    l = 0;
    for j to #A-1 do(
	if (A_j_i == 0) then (
	    sig = append(sig, curList_l);
	    l = l+1;
	    continue;
	    ); 
	sig = append(sig, zeroes);
	);
    sig
    )
