load "DegenMap.m2"
FaceHelperMap = (n1,k1,i) ->(
    A1 = ABuilder(n1,k1);
    A2 = ABuilder(n1-1,k1);
    curList = listBuilder(#A1);
--    zeroes = new List from (binomial(n-1,k)) : 0;
    sig = {};
    zeroes = new List from (#A1) : 0;	 
    remains = {};
    l = 0;
    sig2 = {};
    for j to #A1-1 do(
	if (A1_j_i == 0) then (
	    sig = append(sig, curList_j);
	    continue
	    ); 
	if (A1_j_i == 1) then (
	   remains = append(remains, curList_j); 
	    );
	);
    for j to #sig-1 do(
	 if (A2_j_(i-1) > 0 and l < #remains) then (
	    sig2 = append(sig2, remains_l);
	    l = l+1;
	    continue
	    ); 
    	sig2 = append(sig2,zeroes);
	);
    sig + sig2
    ) 
