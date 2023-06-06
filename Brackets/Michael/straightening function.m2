restart
debug needsPackage "Brackets"

----------------------------------------------------------------straightening function
toBracketPolynomial = method();
toBracketPolynomial(RingElement, BracketRing) := (f, G) -> ( --input: polynomial, bracketring

I = G#ideal;
(f % I) _ G

)

