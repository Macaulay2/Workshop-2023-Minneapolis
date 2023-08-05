needsPackage "CompleteIntersectionResolutions"

Matrix Array := Matrix => (M,L) -> ((M^[L_0])_[L_1])

higherHomotopyFactorization = method();
higherHomotopyFactorization(List,Complex) :=  (L,C) -> (
    Q := ring L_0;
    S := Q[t_1..t_(length L)];
    f := sum(1..length L,i->L_(i-1)*t_i);
    C = C**S;
    H := makeHomotopies(matrix{{f}},chainComplex (C));
    Ln := apply(keys H,i->(i_1+2*(i_0)_0-1,i_1));
    Hn := new HashTable from for i from 0 to length Ln-1 list Ln_i => H#((keys H)#i);
    (lo,hi) := concentration C;
    Lo = select(toList(lo-1..hi+1),i->odd(i));
    Le = select(toList(lo-1..hi+1),i->even(i));
    T1 := table(Lo,Le,(u,v) ->
	if Hn#?(u,v) then map(C_u,C_v,Hn#(u,v)) else 0);
    T2 := table(Le,Lo,(u,v) ->
	if Hn#?(u,v) then map(C_u,C_v,Hn#(u,v)) else 0);
    Co := directSum apply(Lo,i->C_i);
    Ce := directSum apply(Le,i->C_i);
    M1 := map(Ce,Co,matrix T2);
    M2 := map(Co,Ce,matrix T1);
    ZZdfactorization {M1,M2}
    )

higherHomotopyFactorization(RingElement,Complex) := (f,C) -> (
    H := makeHomotopies(matrix{{f}},chainComplex (C));
    Ln := apply(keys H,i->(i_1+2*(i_0)_0-1,i_1));
    Hn := new HashTable from for i from 0 to length Ln-1 list Ln_i => H#((keys H)#i);
    (lo,hi) := concentration C;
    Lo = select(toList(lo-1..hi+1),i->odd(i));
    Le = select(toList(lo-1..hi+1),i->even(i));
    T1 := table(Lo,Le,(u,v) ->
	if Hn#?(u,v) then map(C_u,C_v,Hn#(u,v)) else 0);
    T2 := table(Le,Lo,(u,v) ->
	if Hn#?(u,v) then map(C_u,C_v,Hn#(u,v)) else 0);
    Co := directSum apply(Lo,i->C_i);
    Ce := directSum apply(Le,i->C_i);
    M1 := map(Ce,Co,matrix T2);
    M2 := map(Co,Ce,matrix T1);
    ZZdfactorization {M1,M2}
    )
    
--do higher homotopy factorization coming from a non-resolution
 --do example that contradicts the h(F**F) <= h(F) beta(F)  
    
toBranchedCover = method();
toBranchedCover(ZZdFactorization,Symbol) := (C,z) -> (Q := ring C;
    --code assumes the input is a well-defined facorization
    if not((unique flatten degrees Q)=={0}) then error "Variables from ambient ring should have degree 0";
    d := period C;
    if d==2 then Q.rootOfUnity := -1;
    if Q.?rootOfUnity then t:=Q.rootOfUnity
    else  error "Need to adjoint dth root of unity";
    P := product(d,i->C.dd_i);
    f := sub(P_(0,0),Q);
    S := Q[z];
    zn := (S_*)_0;
    Sk := S/(zn^d+sub(f,S));
    use Sk;
    T := table(toList(0..d-1),toList(0..d-1),(u,v) -> 
	    if u==v then sub(t^(u),Sk)*zn*id_(sub(C_u,Sk)) 
	    else if u==(v-1)%d then sub(C.dd_u,Q)
	    else 0
	    );
    Cterms := directSum for i to d-1 list C_i**Sk;
    map(Cterms,Cterms,sub(matrix T,Sk))
    )

toBranchedCover(ZZdFactorization,RingElement) := (C,z) -> (toBranchedCover(C,getSymbol "z"))


mooreMF = method();
mooreMF(ZZ) := p -> (if p==0 then Q := QQ[a_0..a_2,x_0..x_2];
    if p > 0 then Q := ZZ/p[a_0..a_2,x_0..x_2];
    M1 := matrix{{a_0*x_0,a_1*x_2,a_2*x_1},{a_1*x_1,a_2*x_0,a_0*x_2},{a_2*x_2,a_0*x_1,a_1*x_0}};
    M2 := matrix{{a_1*a_2*x_0^2-a_0^2*x_1*x_2,a_0*a_2*x_1^2-a_1^2*x_0*x_2,a_0*a_1*x_2^2-a_2^2*x_0*x_1},
	{a_0*a_2*x_2^2-a_1^2*x_0*x_1,a_0*a_1*x_0^2-a_2^2*x_1*x_2,a_1*a_2*x_1^2-a_0^2*x_0*x_2},
	{a_0*a_1*x_1^2-a_2^2*x_0*x_2,a_1*a_2*x_2^2-a_0^2*x_0*x_1,a_0*a_2*x_0^2-a_1^2*x_1*x_2}};
    ZZdfactorization {M1,M2}
    )

rk1MCM2gen = (L,d) -> (a := getSymbol "a";
    b := getSymbol "b";
    if d==0 then Q=QQ[x_1..x_4,a,b];
    if d>0 then Q=ZZ/d[x_1..x_4,a,b]; 
    an := (Q_*)_4;
    bn := (Q_*)_5;
    S := Q/(an^2-an+1,bn^2-bn+1);
    use S;
    (i,j,s) := toSequence L;
    M1 := matrix{{x_1-an*x_s,-(x_i^2+bn*x_i*x_j+bn^2*x_j^2)},{x_i-bn*x_j,x_1^2+an*x_1*x_s+an^2*x_s^2}};
    M2 := matrix{{x_1^2+an*x_1*x_s+an^2*x_s^2,x_i^2+bn*x_i*x_j+bn^2*x_j^2},{-(x_i-bn*x_j),x_1-an*x_s}};
    ZZdfactorization {M1,M2}
    )
    	
--need to fix this code	    
rk1MCM3gen = (p,type) -> (a := getSymbol "a";
    b := getSymbol "b";
    c := getSymbol "c";
    d := getSymbol "d";
    e := getSymbol "e";
    if p==0 then Q=QQ[x_1..x_4,a,b,c,d,e,Degrees => {1,1,1,1,0,0,0,0,0}];
    if p>0 then Q=ZZ/p[x_1..x_4,a,b,c,d,e, Degrees => {1,1,1,1,0,0,0,0,0}]; 
    an := (Q_*)_4;
    bn := (Q_*)_5;
    cn := (Q_*)_6;
    dn := (Q_*)_7;
    en := (Q_*)_8;
    S := Q/(an^2-an+1,bn^2-bn+1,cn^2-cn+1,dn^2-dn+1,en^2+en+1,bn*cn*dn-en*an);
    use S;
    if type==1 then (
	M1 := matrix{{0,x_1-an*x_4,x_2-bn*x_3},
	             {x_1-cn*x_2,-bn^2*x_3-an*bn*cn^2*en^2*x_4,bn^2*cn^2*x_3-an*bn*cn*en^2*x_4},
		     {x_3-dn*x_4,cn^2*x_2+bn*cn^2*x_3+an*cn*x_4,-x_1-cn*x_2-an*x_4}};
        M2 := transpose M1;
	return ZZdfactorization {M1,M2};
	);
    if type==2 then (
	M1 := matrix{{0,x_1+x_2,x_3-an*x_4},
	             {x_1+en*x_2,-x_3+cn*x_4,0},
		     {x_3-bn*x_4,0,-x_1-en^2*x_2}};
	M2 := matrix{{0,x_1+x_3,x_2-an*x_4},
	             {x_1-an^2*bn*x_3,-x_2+cn*x_4,0},
		     {x_2-bn*x_4,0,-x_1+an*bn^2*x_3}};
        return ZZdfactorization {M1,M2};
	);
    )

classicalAdjoint = (G) -> (
n := rank target G;
m := rank source G;
matrix table(n, n, (i, j) -> (-1)^(i+j) * det(
submatrix(G, {0..j-1, j+1..n-1},
{0..i-1, i+1..m-1}))));





