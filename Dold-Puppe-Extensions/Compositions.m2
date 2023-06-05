--compositions of n into k parts
posComps = (n,k) -> (sort select(compositions(k,n),i->all(i,j->j>0)))

compHash = (n,k) -> (new HashTable from for i in posComps(n,k)
