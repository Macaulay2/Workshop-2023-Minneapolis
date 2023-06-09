compHash = (n,k) -> (new HashTable from for i in (P=compositions(k,n-k)) list i =>position(P,j->j==i))
    -- Deprecated?

entryCalculator = (mu,j) -> (
    -- Calculates the (mu,j)'th entry of the A matrix
    -- mu is a composition and j ranges from 0 to n
    -- Note the we get compositions with strictly positive parts
    -- by taking compositions with non-negative parts and adding one
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

zeroMatrix = (numRow,numCol) -> toList (numRow:(toList(numCol:0)));
    -- Faster than using table

    
ABuilder = (n,k) -> table(sort compositions(k+1,n-k),0..n,entryCalculator)
    -- This builds the AMatrix
    -- Note that compositions are allowed to have 0 parts (we'll add 1 later)
    -- there are binomial(n,k)-many rows and n+1 columns

rowOfId = (n,l) -> splice {l:0,1,(n-1-l):0}    
-- n-1 so that the result has n-many columns
    

faceMapZero = (n,len) -> (
    -- Output is topDeg-many rows and sum(len,i->binomial(n,i))-many columns
    -* This probably can be optimized
    The output is extremely sparse, so should probably be reimplemented
    as a hash table or function of some kind that specifies rules for 
    the i,j'th entry of the table
    *-
    -- If bounded => true then do this
    -- TODO implement option toggle for bounded complexes
    -- (otherwise why sum these binomial coefficients)
    maxRows = sum(len+1,i->binomial(n-1,i));
    -- The number of rows, taking into account complex length bound
    -- When len is larger than n, we'll have 2^n many rows
    offset = 0; -- keeps track of what the top zero pad should be
    verticalStrips = for k to maxK list (
    -- First loop over the vertical strips of columns
	bink=binomial(n,k); -- Each strip is width binomial(n,k)
    	topZeros = zeroMatrix(offset,bink);
    -- Create appropriately sized top pad of zeros
	modifiedMat = for row to bink-1 list (
    -- This modified identity matrix uses Pascal's identity
    -- First binomial(n-1,k-1) entries will be differential
    	    if row<=binomial(n-1,k-1)-1 then splice{row:0,k+1,(bink-1-row):0}
    --k+1 to keep track of which differential we'll need (k'th diff of complex)
	    else splice{row:0, -k-1,(bink-1-row):0}
    -- -k-1 to keep track of which identity map we'll need (identity map of C_k)
	    );
	botZeros = zeroMatrix(maxRows-offset-binomial(n,k),binomial(n,k));
    -- Create appropriately sized bottom pad of zeros
    	offset = offset + binomial(n-1,k-1);
    -- Update the offset
    	verticalStrip = flatten{topZeros,modifiedMat,botZeros}
    -- Concat the modified matrix with the pads
	);
    for row to maxRows-1 list (
    -- This concats the rows of the vertical strips together
    	flatten apply(verticalStrips, strip->strip#row)
	)
    )


faceMapnk = (n,k) -> (
    -- output is binomial(n-1,k)-many rows and binomial(n,k)-many columns
    if k==n then 
    -- I don't remember if we actually need this check
    -- I think we need to zero pad extra rows at when faceMapi is called below
	{{0}}
    else (
	Abigvector = apply(ABuilder(n,k), row->row_n);
    -- Selects the nth column of AMat(n,k)
	Asmallvector = apply(ABuilder(n-1,k), row->row_(n-1));
    -- Selects the (n-1)th column of AMat(n-1,k)
	onePos = positions(Abigvector, i -> i==1);
    -- The positions of ones in Asmallvector describe horizontal locations where we'll want a row of the identity matrix
	outputMat = onePos / (l -> rowOfId(binomial(n,k),l))
	)
    )
 
faceMapik = (n,k,i) -> (
    -- I think there are supposed to be binomial(n-1,k)-many rows
    -- There are binomial(n,k)-many columns
    -- I think my timing code was wack and this is actually slower than the earlier implementation
    -- Need to go back and check that
    -- From what I recall, they were both pretty close though so it's not a high priority
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

faceMapi = (n,i,C) -> (
    -- Need to check if this can be sped up by direct summing as lists
    -- Instead of first passing to matrices
    -- and then saving the creation of the matrix option till the very end
    maxK = min(length C,n-1);
    -- Are we running till maxK? or until maxK+1?
    if i==0 then (
	maxK = min(length C,n);
	promoteFaceMapZerotoComplex(faceMapZero(n,maxK),C)
	)
    else if i==n then (
	preMat = fold(directSum,for k from 0 to maxK list (
	    promoteMaptoComplex(faceMapnk(n,k),k,C)
    	    -- Notice that we only pass in k at most n-1, so the check for k==n in faceMapi isn't needed
	    )); 
	preMat | map(target preMat,C_n,0) 
	)
    else (
    	preMat = fold(directSum,for k from 0 to maxK list (
	     promoteMaptoComplex(faceMapik(n,k,i),k,C)
	    )); 
	preMat | map(target preMat,C_n,0) 
	)
    )

degenMapik = (n, k, i) -> (
    -- We'll have binomial(n,k) many columns #A(n+1,k)=binomial(n+1,k)-many rows
    
    numCols = binomial(n,k);
    zeroes = new List from (numCols) : 0;
    A = ABuilder(n+1,k);
    sig = new MutableList from 0..(#A-1);
    l = 0;
    for j to #A-1 do(
	if (A_j_i == 0) then (
	    sig#j = rowOfId(numCols,l);
	    l = l+1;
	    continue;
	    ); 
	sig#j = zeroes;
	);
    new List from sig
    )

degenMapi = (n,i,C) -> (
    maxK = min(length C,n);
    preMat = fold(directSum,for k from 0 to maxK list (
	    promoteMaptoComplex(degenMapik(n,k,i),k,C)
	    )
	);
    preMat || map(C_(n+1),source preMat,0)
    )
    
promoteMaptoComplex = (d,k,C) -> (matrix d)**id_(C_k)
promoteFaceMapZerotoComplex = (d,C) -> matrix (
   -- This runs over the whole matrix again and so is probably not an efficient way of doing things
   -- Having the output of FaceMapZero be a hash table would be faster
   d / (i -> apply(i,j -> (
	       if j<0 then (
	           id_(C_(abs(j+1)))
    	    	    -- abs because identities were stored as negatives
		    -- +1 first because we stored identity of C_k as -k-1
	       ) 
	       else if j>0 then (
		   dd^C_(j-1)
    	    	    -- -1 because we stored the differential of C_k as k+1
		    -- I don't actually think we need to do that though?  probably could just store it as k
		   ) 
	       else 0
	    )
	)
    )
   )
	   
	   
	   
	   
	   
	   
	   
	    
