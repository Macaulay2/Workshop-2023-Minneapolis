doc ///
    Key
        (ASMFullList, ZZ)
        ASMFullList
    Headline
        lists all ASMs of a fixed size
    Usage
        ASMFullList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $1 \leq n \leq 7$, outputs the full list of $n \times n$ alternating sign matrices,
	    including permutation matrices.
///

doc ///
    Key
        (ASMRandomList, ZZ,ZZ)
        ASMRandomList
    Headline
        lists random ASMs of a fixed size
    Usage
        ASMFullList(n,m)
    Inputs
        n:ZZ
        m:ZZ
    Description
        Text
            For $1 \leq n \leq 7$, this function lists $m$ random $n\times n$ alternating sign matrices.
///

doc ///
    Key
        (cohenMacaulayASMsList, ZZ)
        cohenMacaulayASMsList
    Headline
        lists all Cohen-Macaulay ASMs of a fixed size which are not permutation matrices
    Usage
        cohenMacaulayASMsList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $1 \leq n \leq 6$, this function lists all $n\times n$ alternating sign matrices $A$ which are \textbf{not} permutation matrices
	    such that the corresponding ASM variety is Cohen-Macaulay.
	    By a theorem of Fulton [Ful92], permutation matrices always have Cohen-Macaulay Schubert determinantal ideals.
	    
	    @UL {
	    {"[Ful92] William Fulton, ",
	    HREF("https://sites.math.washington.edu/~billey/classes/schubert.library/fulton.essential.set.pdf",
		EM "Flags, Schubert polynomials, degeneracy loci, and determinantal formulas"),
	    " , Duke Math J. 65 (1992): 381-420."}
	    }@
///

doc ///
    Key
        (nonCohenMacaulayASMsList, ZZ)
        nonCohenMacaulayASMsList
    Headline
        lists all non-Cohen-Macaulay ASMs of a fixed size
    Usage
        nonCohenMacaulayASMsList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $1 \leq n \leq 6$, this function lists all ASMs of size $n$ such that the corresponding ASM variety
	    is \textbf{not} Cohen-Macaulay.
///

doc ///
    Key
        (initialIdealsList, ZZ)
        initialIdealsList
    Headline
        lists all antidiagonal initial ideals of ASMs of a fixed size
    Usage
        initialIdealsList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $3 \leq n \leq 6$, this function lists all antidiagonal initial ideals for ASMs of size $n$.
///
