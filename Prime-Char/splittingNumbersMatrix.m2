expDecomp = (R,p,e,f) -> apply(exponents(f),exponent->{coefficient(R_exponent,f)*R_(exponent //p^e),exponent%p^e});
--Gets the exponent vectors of each monomial X^u of the polynomial f, and associates to u the two-element list whose
        --first entry is cX^v and second entry is w, where c is the coefficient of X^u in f and u = p^e*v + w.
        

AKmatrix = (R,e,f) -> (
    p:=char R;
    d:=dim R;
    n:=p^e;
    I:=(ideal(vars R))^[n];
    monBasis:=(entries basis(R^1/I))#0;
    expBasis:=flatten(apply(monBasis,exponents));
    L:=for i from 0 to p^(e*d)-1 list f*(monBasis#i);
    T:=(for i from 0 to p^(e*d)-1 list expDecomp(R,p,e,L#i));
    return transpose matrix(for t in T list( for exponent in expBasis list sum(
	apply(select(t, element -> (element#1 == exponent)),j->j#0))))
)


AKmatrixIdeal = (R,e,I) -> (
    p:=char R;
    d:=dim R;
    n:=p^e;
    M:=fold((i,j) -> (i|j), for g in I_* list AKmatrix(R,e,g));
    return M
    )
    
    
