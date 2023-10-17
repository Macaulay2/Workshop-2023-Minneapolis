doc ///
    Key
        (ASMFullList, ZZ)
        ASMFullList
    Headline
        provides a list of all ASMs of a certain size
    Usage
        ASMFullList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $1 \leq n \leq 7$, a full list of all ASMs of size $n$ comes with the package.
///

doc ///
    Key
        (ASMRandomList, ZZ,ZZ)
        ASMRandomList
    Headline
        provides a list of some random ASMs of a certain size
    Usage
        ASMFullList(n,m)
    Inputs
        n:ZZ
        m:ZZ
    Description
        Text
            For $1 \leq n \leq 7$, a list of $m$ random ASMs of size $n$ is outputted.
///

doc ///
    Key
        (cohenMacaulayASMsList, ZZ)
        cohenMacaulayASMsList
    Headline
        provides a list of all Cohen-Macaulay ASMs which are not permutation matrices of a certain size
    Usage
        cohenMacaulayASMsList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $1 \leq n \leq 6$, a list of all Cohen-Macaulay ASMs which are not permutation matrices of size $n$ comes with the package.
	    It is known by Fulton [Ful92] that the permutation matrices have Cohen-Macaulay Schubert determinantal ideals.
	    
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
        provides a list of all non-Cohen-Macaulay ASMs of a certain size
    Usage
        nonCohenMacaulayASMsList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $1 \leq n \leq 6$, a list of all non-Cohen-Macaulay ASMs of size $n$ comes with the package.
///

doc ///
    Key
        (initialIdealsList, ZZ)
        initialIdealsList
    Headline
        provides a list of all initial ideals of ASMs of a certain size
    Usage
        initialIdealsList(n)
    Inputs
        n:ZZ
    Description
        Text
            For $3 \leq n \leq 6$, a list of all initial ideals of ASMs of size $n$ comes with the package.
///
