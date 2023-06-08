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

zeroMatrix = (numRow,numCol) -> toList (numRow:(toList(numCol:0)));

    
ABuilder = (n,k) -> table(sort compositions(k+1,n-k),0..n,entryCalculator)

rowOfId = (n,l) -> splice {l:0,1,(n-1-l):0}    

    

faceMapZero = (n,len) -> (
    -* This probably can be optimized
    The output is extremely sparse, so should probably be reimplemented
    as a hash table or function of some kind that specifies rules for 
    the i,j'th entry of the table
    *-
    -- If bounded => true then do this
    -- TODO implement option toggle for bounded complexes
    topDeg = sum(len+1,i->binomial(n-1,i));
    offset = 0; -- keeps track of where we are
    verticalStrips = for k to len list (
    -- First loop over the vertical strips of columns
    -- Create appropriately sized top pad of zeros
	bink=binomial(n,k);
    	topZeros = zeroMatrix(offset,bink);
	modifiedMat = for row to bink-1 list (
    	    if row<=binomial(n-1,k-1)-1 then splice{row:0,k+1,(bink-1-row):0}
    --k+1 to keep track of which differential we'll need
	    else splice{row:0, -k-1,(bink-1-row):0}
    -- -k-1 to keep track of which identity map we'll need
	    );
    -- Create appropriately sized bottom pad of zeros
	botZeros = zeroMatrix(topDeg-offset-binomial(n,k),binomial(n,k));
    	offset = offset + binomial(n-1,k-1);
    	verticalStrip = flatten{topZeros,modifiedMat,botZeros}
    -- Update Offset
	);
    for row to topDeg-1 list (
    	flatten apply(verticalStrips, strip->strip#row)
	)
    )


faceMapnk = (n,k) -> (
    if k==n then 
	{{0}} 
    else
	Abigvector = apply(ABuilder(n,k), row->row_n);
	Asmallvector = apply(ABuilder(n-1,k), row->row_(n-1));
	onePos = positions(Abigvector, i -> i==1);
	outputMat = onePos / (l -> rowOfId(binomial(n,k),l))
    )
 
faceMapik = (n,k,i) -> (
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
    len = length C;
    -- Need to check if this can be sped up by direct summing as lists
    -- Instead of first passing to matrices
    -- and then saving the creation of the matrix option till the very end
    if i==0 then promoteFaceMapZerotoComplex(faceMapZero(n,len),C)
    else if i==n then (
	fold(directSum,for k from 0 to len list (
	    promoteMaptoComplex(faceMapnk(n,k),k,C)
	    ))
	)
    else (
    	fold(directSum,for k from 0 to len list (
	     promoteMaptoComplex(faceMapik(n,k,i),k,C)
	     ))
	)
    )

degenMapik = (n, k, i) -> (
    idSize = binomial(n-2,k);
    zeroes = new List from (binomial(n-2,k)) : 0;
    A = ABuilder(n-1,k);
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

degenMapi = (n,i,C) -> (
    len = length C;
    fold(directSum,for k from 0 to len list (
	    promoteMaptoComplex(degenMapik(n,k,i),k,C)
	    )
	)
    )
    
promoteMaptoComplex = (d,k,C) -> (matrix d)**id_(C_k)
promoteFaceMapZerotoComplex = (d,C) -> matrix (
   -- This runs over the whole matrix again and so is probably not an efficient way of doing things
   -- Having the output of FaceMapZero be a hash table would be faster
   d / (i -> apply(i,j -> (
	       if j<0 then (
	           id_(C_(abs(j+1)))
	       ) 
	       else if j>0 then (
		   dd^C_(j-1)
		   ) 
	       else 0
	    )
	)
    )
   )
	   
	   
	   
	   
	   
	   
	   
	    
