newPackage(
    "FInvariants",
    Version => "1.1",
    Date => "June 8, 2023",
    Authors => {{Name => "Anna Brosowsky", Email => "annabro@umich.edu"},
	{Name => "Havi Ellers", Email => "ellers@umich.edu"},
	{Name => "Moty Katzman", Email => "M.Katzman@shef.ac.uk"},
	{Name => "Jiamin Li", Email => "jli283@uic.edu"},
	{Name => "Swaraj Pande", Email => "swarajsp@umich.edu"},
	{Name => "Abraham Pascoe", Email => "abrahampascoe@ku.edu"},
	{Name => "Austyn Simpson", Email => "austyn@umich.edu"},
	{Name => "Pedro Teixeira", Email => "pteixeir@knox.edu"}
        },
    PackageImports => {"TestIdeals"},
    Headline => "A package for calculations of Frobenius invariants.",
    AuxiliaryFiles => true
    )

export {"HSL"}

HSL = method()


HSL (Module,Matrix) := (A,u) -> (
    A1:=matrix(entries(presentation(A)));
    u1:=matrix(entries(u));
    e:=0;
    M:=gens(target(u1));
    while true do (
	newM:=matrix(entries(frobeniusRoot(1,u1*M)));
	if image(newM)+image(A1) == image(M)+image(A1) then return e;
	e=e+1;
	M=newM;
	)
    )


HSL (Matrix,Matrix) := (A,u) -> (
    HSL(coker(A),u)
    )

beginDocumentation()



end--


--Test example, should output 2

restart
needsPackage "FInvariants"
needsPackage "TestIdeals"
R=ZZ/2[x1,x2,x3,x4,x5]
I=ideal((x2)^2+(x1)*(x3),(x1)*(x2)*(x4)^2+(x3)^3*(x5),(x1)^2*(x4)^2+(x2)*(x3)^2*(x5))
M2=R^1/I
A=Ext^2(M2,R^1)
M1=R^1/frobenius(1,I)
phi=map(M2,M1,1)
u=Ext^2(phi,R^1)
B=presentation(A)
HSL(B,u)


--Test example, should output 1

restart
needsPackage "FInvariants"
needsPackage "TestIdeals"
R=ZZ/11[x,y,z]
I=ideal(x^5+y^5+z^5)
M2=R^1/I
Ext(M2,R^1)
A=Ext^1(M2,R^1)
M1=R^1/frobenius(1,I)
phi=map(M2,M1,1)
u=Ext^1(phi,R^1)
B=presentation(A)
HSL(B,u)
